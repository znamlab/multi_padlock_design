"""NUPACK-based exhaustive heterodimer analysis for padlock probes.

This module enumerates all unique unordered probe (padlock) pairs, batches the
pairs across SLURM jobs (via the slurm_it decorator), and for each pair computes:
  * ΔG (kcal/mol) of the heterodimer complex (from NUPACK partition function)
  * Equilibrium concentration of the heterodimer complex
  * Fraction / percent of each monomer present in heterodimer form relative to
    its initial molar concentration (padlock_conc).

Two matrices are produced by the aggregation step:
  1. dg_matrix (NxN): symmetric matrix of heterodimer free energies (kcal/mol)
  2. percent_matrix (NxN): symmetric matrix of percent heterodimer bound

Files written:
  heterodimer_nupack_pairs_part_XXXX.pkl   (batch part files)
  heterodimer_nupack_dg_matrix.pkl         (pickle of {matrix, names, ...})
  heterodimer_nupack_percent_matrix.pkl    (pickle of {matrix, names, ...})
  heterodimer_nupack_results.pkl           (single pickle bundling both matrices)

Example (submitted from a notebook or script):

    from lib import nupack_heterodimer
    nupack_heterodimer.submit_exhaustive_heterodimer_jobs(
        probe_df_path="monahan_panel_barcoded.csv",
        output_dir="/camp/lab/znamenskiyp/scratch/olfr_dg_hetero",
        slurm_folder=str(Path.home()/"slurm_logs"/"olfr_probe_design"),
        n_batches=900,
        time="3-00:00:00",
        mem="8G",
        partition="ncpu",
        sequence_col="padlock",
        dependency_aggregate_time="02:00:00",
        dependency_aggregate_mem="32G",
        dry_run=False,
    )

Dependencies: nupack, pandas, numpy, tqdm (optional for progress bars).
"""

from __future__ import annotations

import importlib
import math
import os  # added for sequence/name caching
import pickle
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from multi_padlock_design.secondary_structure.check_padlocks import clean_seq, slurm_it

try:  # NUPACK imports
    _nupack = importlib.import_module("nupack")
except ImportError as e:
    raise ImportError(
        "NUPACK is required for nupack_heterodimer module. Install nupack first."
    ) from e

Model = _nupack.Model
SetSpec = _nupack.SetSpec
Strand = _nupack.Strand
Tube = _nupack.Tube
tube_analysis = _nupack.tube_analysis

# ------------------------------- Pair Enumeration -------------------------------


def pair_generator(N: int, start_k: int, end_k: int) -> Iterable[Tuple[int, int]]:
    """Yield (i, j) index pairs for linearized upper-triangular indices.

    Pairs are ordered lexicographically (i, j) for 0 <= i < j < N.
    start_k and end_k bound the half-open interval of the linear index space to
    enumerate for the current batch.
    """
    remaining = start_k
    i = 0
    while i < N - 1 and remaining >= (N - i - 1):
        remaining -= N - i - 1
        i += 1
    if i >= N - 1:
        return
    j = i + 1 + remaining
    k = start_k
    while k < end_k and i < N - 1:
        yield i, j
        k += 1
        j += 1
        if j >= N:
            i += 1
            j = i + 1


# --------------------------- Core Heterodimer Computation ---------------------------


def _compute_pair_nupack(
    s1: str,
    s2: str,
    name1: str,
    name2: str,
    padlock_conc: float,
    model_kwargs: Dict[str, Any],
    model: Optional[Model] = None,
) -> Dict[str, Any]:
    """Compute heterodimer ΔG and fraction bound for a single pair.

    Returns dict with keys: heterodimer_dG, heterodimer_M, heterodimer_fraction,
    heterodimer_percent. On failure, values are math.nan.
    """
    import math as _math

    if model is None:
        model = Model(**model_kwargs)

    strand1 = Strand(s1, name=name1)
    strand2 = Strand(s2, name=name2)
    tube = Tube(
        strands={strand1: padlock_conc, strand2: padlock_conc},
        complexes=SetSpec(max_size=2),
        name=f"Dimer_{name1}_{name2}",
    )
    res_map = tube_analysis(tubes=[tube], compute=["pfunc"], model=model)
    tube_res = res_map[tube]
    complex_concs = getattr(tube_res, "complex_concentrations", {})
    heterodimer_conc = _math.nan
    for cx, conc in complex_concs.items():  # locate 2-strand complex
        try:
            names = [s.name for s in cx.strands]
        except AttributeError:  # fallback if structure differs
            names = [s.name for s in cx]
        if len(names) == 2 and name1 in names and name2 in names:
            heterodimer_conc = conc
            break
    dG = _math.nan
    for key in (f"({name1}+{name2})", f"({name2}+{name1})"):
        if key in res_map:
            dG = getattr(res_map[key], "free_energy", _math.nan)
            break
    if padlock_conc > 0 and not _math.isnan(heterodimer_conc):
        fraction = heterodimer_conc / padlock_conc
        percent = fraction * 100.0
    else:
        fraction = _math.nan
        percent = _math.nan
    return {
        "heterodimer_dG": dG,
        "heterodimer_M": heterodimer_conc,
        "heterodimer_fraction": fraction,
        "heterodimer_percent": percent,
    }


# --------------------------------- Batch Function ---------------------------------


@slurm_it(
    conda_env="multi-padlock-design",
    from_imports={"lib.nupack_heterodimer": "run_nupack_heterodimer_batch"},
)
def run_nupack_heterodimer_batch(
    probe_df_path: str,
    output_dir: str,
    batch_index: int,
    n_batches: int = 1000,
    sequence_col: str = "padlock",
    name_col: Optional[str] = None,
    padlock_conc: float = 1e-9,
    overwrite: bool = False,
    use_seq_cache: bool = True,
    seq_cache_name: str = "probe_sequences_cache.pkl",
) -> str:
    """Compute NUPACK heterodimer metrics for one batch of unordered probe pairs.

    Adds optional on-disk caching of cleaned sequences & names so that many
    concurrent batch jobs avoid repeatedly reading / parsing a large CSV.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"heterodimer_nupack_pairs_part_{batch_index:04d}.pkl"
    if out_file.exists() and not overwrite:
        return str(out_file)

    probe_df_path = Path(probe_df_path)
    cache_file = out_dir / seq_cache_name
    loaded_from_cache = False
    if use_seq_cache and cache_file.exists():
        try:
            with cache_file.open("rb") as f:
                cache_obj = pickle.load(f)
            seqs = cache_obj["seqs"]
            names = cache_obj["names"]
            loaded_from_cache = True
        except Exception:
            loaded_from_cache = False

    if not loaded_from_cache:
        if probe_df_path.suffix.lower() in {".pkl", ".pickle"}:
            df = pd.read_pickle(probe_df_path)
        else:
            df = pd.read_csv(probe_df_path)
        if sequence_col not in df.columns:
            raise ValueError(f"Sequence column '{sequence_col}' not found")
        if name_col is None:
            for cand in ["name", "probe_name", "id", "gene", "padlock_name"]:
                if cand in df.columns:
                    name_col = cand
                    break
        if name_col is None:
            name_col = "_autoname"
            df = df.copy()
            df[name_col] = [f"probe_{i}" for i in range(len(df))]
        seqs = df[sequence_col].map(clean_seq).tolist()
        names = df[name_col].astype(str).tolist()
        if use_seq_cache:
            tmp_file = cache_file.with_suffix(".tmp")
            try:
                with tmp_file.open("wb") as f:
                    pickle.dump(
                        {"seqs": seqs, "names": names, "sequence_col": sequence_col}, f
                    )
                os.replace(tmp_file, cache_file)
            except Exception:
                if tmp_file.exists():
                    try:
                        tmp_file.unlink()
                    except Exception:
                        pass

    N = len(seqs)
    if N < 2:
        raise ValueError("Need at least two probes")
    total_pairs = N * (N - 1) // 2
    n_batches = min(n_batches, total_pairs) if total_pairs else 0
    if n_batches == 0:
        raise ValueError("No pairs to enumerate")
    if batch_index >= n_batches:
        with out_file.open("wb") as f:
            pickle.dump(pd.DataFrame([]), f)
        return str(out_file)

    per_batch = math.ceil(total_pairs / n_batches)
    start_k = batch_index * per_batch
    end_k = min(start_k + per_batch, total_pairs)
    if start_k >= end_k:
        with out_file.open("wb") as f:
            pickle.dump(pd.DataFrame([]), f)
        return str(out_file)

    model_kwargs = dict(material="dna", celsius=45, sodium=0.075, magnesium=0.01)
    model = Model(**model_kwargs)  # reuse same model for all pairs in this batch

    from tqdm.auto import tqdm

    batch_total = end_k - start_k
    miniters = max(1, batch_total // 200) if batch_total > 0 else 1
    pbar = tqdm(
        total=batch_total,
        desc=f"nupack heterodimer batch {batch_index+1}/{n_batches}",
        unit="pair",
        smoothing=0.0,
        mininterval=30.0,
        miniters=miniters,
        dynamic_ncols=False,
        disable=False,
        ascii=True,
        leave=False,
    )

    records: List[Dict[str, Any]] = []
    for i, j in pair_generator(N, start_k, end_k):
        s1 = seqs[i]
        s2 = seqs[j]
        name1 = names[i]
        name2 = names[j]
        metrics = _compute_pair_nupack(
            s1,
            s2,
            name1,
            name2,
            padlock_conc,
            model_kwargs,
            model=model,
        )
        records.append({"i": i, "j": j, "probe_i": name1, "probe_j": name2, **metrics})
        pbar.update(1)
    pbar.close()

    part_df = pd.DataFrame(records)
    with out_file.open("wb") as f:
        pickle.dump(part_df, f)
    return str(out_file)


# -------------------------------- Aggregation Function --------------------------------


@slurm_it(
    conda_env="multi-padlock-design",
    from_imports={"lib.nupack_heterodimer": "aggregate_nupack_heterodimer_results"},
)
def aggregate_nupack_heterodimer_results(
    output_dir: str,
    probe_df_path: str,
    sequence_col: str = "padlock",
    dg_matrix_pickle: Optional[str] = None,
    percent_matrix_pickle: Optional[str] = None,
    combined_pickle: Optional[str] = None,
    overwrite: bool = True,
    use_seq_cache: bool = True,
    seq_cache_name: str = "probe_sequences_cache.pkl",
) -> str:
    """Aggregate part files into symmetric ΔG and percent matrices.

    Will reuse sequence/name cache if available to avoid large CSV re-reads.
    """
    out_dir = Path(output_dir)
    if dg_matrix_pickle is None:
        dg_matrix_pickle = out_dir / "heterodimer_nupack_dg_matrix.pkl"
    else:
        dg_matrix_pickle = Path(dg_matrix_pickle)
    if percent_matrix_pickle is None:
        percent_matrix_pickle = out_dir / "heterodimer_nupack_percent_matrix.pkl"
    else:
        percent_matrix_pickle = Path(percent_matrix_pickle)
    if combined_pickle is None:
        combined_pickle = out_dir / "heterodimer_nupack_results.pkl"
    else:
        combined_pickle = Path(combined_pickle)

    if (
        dg_matrix_pickle.exists()
        and percent_matrix_pickle.exists()
        and combined_pickle.exists()
        and not overwrite
    ):
        return str(combined_pickle)

    probe_df_path = Path(probe_df_path)
    cache_file = out_dir / seq_cache_name
    names: Optional[List[str]] = None
    if use_seq_cache and cache_file.exists():
        try:
            with cache_file.open("rb") as f:
                cache_obj = pickle.load(f)
            if cache_obj.get("sequence_col") == sequence_col:
                names = cache_obj["names"]
        except Exception:
            names = None

    if names is None:  # fallback to reading original source
        if probe_df_path.suffix.lower() in {".pkl", ".pickle"}:
            df = pd.read_pickle(probe_df_path)
        else:
            df = pd.read_csv(probe_df_path)
        if sequence_col not in df.columns:
            raise ValueError(f"Sequence column '{sequence_col}' not found")
        name_col = None
        for c in ["name", "probe_name", "id", "gene", "padlock_name", "_autoname"]:
            if c in df.columns:
                name_col = c
                break
        if name_col is None:
            name_col = "_autoname"
            df[name_col] = [f"probe_{i}" for i in range(len(df))]
        names = df[name_col].astype(str).tolist()

    N = len(names)
    dg_matrix = np.full((N, N), np.nan, dtype=float)
    percent_matrix = np.full((N, N), np.nan, dtype=float)

    part_files = sorted(out_dir.glob("heterodimer_nupack_pairs_part_*.pkl"))
    if not part_files:
        raise FileNotFoundError(
            "No part files found to aggregate (pattern "
            "heterodimer_nupack_pairs_part_*.pkl)"
        )

    # Count total rows for progress
    total_rows = 0
    for pf in part_files:
        try:
            with pf.open("rb") as f:
                part_df = pickle.load(f)
        except Exception:
            continue
        if isinstance(part_df, pd.DataFrame):
            total_rows += len(part_df)

    try:
        from tqdm.auto import tqdm
    except ImportError:

        def tqdm(x, **k):
            return x

    pbar = tqdm(
        total=total_rows,
        desc="aggregate nupack heterodimer",
        unit="pair",
        mininterval=30.0,
        smoothing=0.0,
        ascii=True,
        disable=False,
        leave=False,
    )

    for pf in part_files:
        try:
            with pf.open("rb") as f:
                part_df = pickle.load(f)
        except Exception:
            continue
        if not isinstance(part_df, pd.DataFrame):
            continue
        for row in part_df.itertuples():
            i = row.i
            j = row.j
            if 0 <= i < N and 0 <= j < N:
                dg = getattr(row, "heterodimer_dG", math.nan)
                pct = getattr(row, "heterodimer_percent", math.nan)
                dg_matrix[i, j] = dg_matrix[j, i] = dg
                percent_matrix[i, j] = percent_matrix[j, i] = pct
            pbar.update(1)
    pbar.close()

    np.fill_diagonal(dg_matrix, 0.0)  # self ΔG set to 0.0 for convenience
    np.fill_diagonal(percent_matrix, 0.0)  # self percent = 0 (no heterodimer)

    with open(dg_matrix_pickle, "wb") as f:
        pickle.dump(
            {"matrix": dg_matrix, "names": names, "sequence_col": sequence_col}, f
        )
    with open(percent_matrix_pickle, "wb") as f:
        pickle.dump(
            {"matrix": percent_matrix, "names": names, "sequence_col": sequence_col}, f
        )
    with open(combined_pickle, "wb") as f:
        pickle.dump(
            {
                "dg_matrix": dg_matrix,
                "percent_matrix": percent_matrix,
                "names": names,
                "sequence_col": sequence_col,
            },
            f,
        )
    return str(combined_pickle)


# --------------------------------- Job Submission ---------------------------------


def submit_exhaustive_heterodimer_jobs(
    probe_df_path: str,
    output_dir: str,
    slurm_folder: str,
    n_batches: int = 1000,
    time: str = "12:00:00",
    mem: str = "8G",
    partition: Optional[str] = "ncpu",
    cpus_per_task: int = 1,
    sequence_col: str = "padlock",
    dependency_aggregate_time: str = "02:00:00",
    dependency_aggregate_mem: str = "32G",
    padlock_conc: float = 1e-9,
    dry_run: bool = False,
    use_seq_cache: bool = True,
    seq_cache_name: str = "probe_sequences_cache.pkl",
) -> tuple[list[str], Optional[str]]:
    """Submit batched NUPACK heterodimer jobs plus a dependent aggregation job.

    Sequence/name caching is enabled by default to reduce concurrent filesystem load.
    """
    probe_df_path_p = Path(probe_df_path)
    if probe_df_path_p.suffix.lower() in {".pkl", ".pickle"}:
        df = pd.read_pickle(probe_df_path_p)
    else:
        df = pd.read_csv(probe_df_path_p)
    if sequence_col not in df.columns:
        raise ValueError(f"Sequence column '{sequence_col}' not found")
    N = len(df)
    if N < 2:
        raise ValueError("Need at least two probes")
    total_pairs = N * (N - 1) // 2
    effective_batches = min(n_batches, total_pairs)
    if effective_batches < n_batches:
        n_batches = effective_batches

    if dry_run:
        print(f"[DRY RUN] N={N} total_pairs={total_pairs} n_batches={n_batches}")
        return [], None

    slurm_opts_batch = {"time": time, "mem": mem, "cpus-per-task": cpus_per_task}
    if partition:
        slurm_opts_batch["partition"] = partition
    slurm_opts_agg = {
        "time": dependency_aggregate_time,
        "mem": dependency_aggregate_mem,
        "cpus-per-task": 1,
    }
    if partition:
        slurm_opts_agg["partition"] = partition

    batch_job_ids: list[str] = []
    for i in range(n_batches):
        jid = run_nupack_heterodimer_batch(
            probe_df_path=str(probe_df_path_p),
            output_dir=str(output_dir),
            batch_index=i,
            n_batches=n_batches,
            sequence_col=sequence_col,
            padlock_conc=padlock_conc,
            slurm_folder=slurm_folder,
            scripts_name=f"nupack_heterodimer_batch_{i:04d}",
            slurm_options=slurm_opts_batch,
            use_slurm=True,
            use_seq_cache=use_seq_cache,
            seq_cache_name=seq_cache_name,
        )
        batch_job_ids.append(jid)

    aggregate_job_id = aggregate_nupack_heterodimer_results(
        output_dir=str(output_dir),
        probe_df_path=str(probe_df_path_p),
        sequence_col=sequence_col,
        slurm_folder=slurm_folder,
        scripts_name="nupack_heterodimer_aggregate",
        job_dependency=batch_job_ids,
        slurm_options=slurm_opts_agg,
        use_slurm=True,
        use_seq_cache=use_seq_cache,
        seq_cache_name=seq_cache_name,
    )
    return batch_job_ids, aggregate_job_id


# -------------------------------- Selected-vs-All (Per-Probe) -------------------------


@slurm_it(
    conda_env="multi-padlock-design",
    from_imports={"lib.nupack_heterodimer": "run_nupack_selected_vs_all_for_probe"},
)
def run_nupack_selected_vs_all_for_probe(
    probe_df_path: str,
    output_dir: str,
    probe_index: int,
    sequence_col: str = "padlock",
    name_col: Optional[str] = None,
    padlock_conc: float = 1e-9,
    overwrite: bool = False,
    use_seq_cache: bool = True,
    seq_cache_name: str = "probe_sequences_cache.pkl",
) -> str:
    """Compute NUPACK heterodimer metrics for one selected probe vs all probes.

    Writes a part file with rows for j in [0..N-1] and columns:
      i, j, probe_i, probe_j, heterodimer_dG, heterodimer_M,
      heterodimer_fraction, heterodimer_percent
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"selected_vs_all_probe_{probe_index:06d}.pkl"
    if out_file.exists() and not overwrite:
        return str(out_file)

    probe_df_path = Path(probe_df_path)
    cache_file = out_dir / seq_cache_name

    # Load cached seqs/names when available
    seqs: List[str]
    names: List[str]
    loaded_from_cache = False
    if use_seq_cache and cache_file.exists():
        try:
            with cache_file.open("rb") as f:
                cache_obj = pickle.load(f)
            seqs = cache_obj["seqs"]
            names = cache_obj["names"]
            loaded_from_cache = True
        except Exception:
            loaded_from_cache = False

    if not loaded_from_cache:
        if probe_df_path.suffix.lower() in {".pkl", ".pickle"}:
            df = pd.read_pickle(probe_df_path)
        else:
            df = pd.read_csv(probe_df_path)
        if sequence_col not in df.columns:
            raise ValueError(f"Sequence column '{sequence_col}' not found")
        if name_col is None:
            for cand in ["name", "probe_name", "id", "gene", "padlock_name"]:
                if cand in df.columns:
                    name_col = cand
                    break
        if name_col is None:
            name_col = "_autoname"
            df = df.copy()
            df[name_col] = [f"probe_{i}" for i in range(len(df))]
        seqs = df[sequence_col].map(clean_seq).tolist()
        names = df[name_col].astype(str).tolist()
        if use_seq_cache:
            tmp_file = cache_file.with_suffix(".tmp")
            try:
                with tmp_file.open("wb") as f:
                    pickle.dump(
                        {"seqs": seqs, "names": names, "sequence_col": sequence_col}, f
                    )
                os.replace(tmp_file, cache_file)
            except Exception:
                if tmp_file.exists():
                    try:
                        tmp_file.unlink()
                    except Exception:
                        pass

    N = len(seqs)
    if not (0 <= probe_index < N):
        raise IndexError(f"probe_index {probe_index} out of range for N={N}")

    model_kwargs = dict(material="dna", celsius=45, sodium=0.075, magnesium=0.01)
    model = Model(**model_kwargs)

    from tqdm.auto import tqdm  # lazy import in job for progress (optional)

    records: List[Dict[str, Any]] = []
    pbar = tqdm(
        total=N,
        desc=f"nupack selected-vs-all i={probe_index}",
        unit="pair",
        ascii=True,
        disable=False,
        leave=False,
    )
    s1 = seqs[probe_index]
    name1 = names[probe_index]
    for j in range(N):
        s2 = seqs[j]
        name2 = names[j]
        metrics = _compute_pair_nupack(
            s1=s1,
            s2=s2,
            name1=name1,
            name2=name2,
            padlock_conc=padlock_conc,
            model_kwargs=model_kwargs,
            model=model,
        )
        records.append(
            {"i": probe_index, "j": j, "probe_i": name1, "probe_j": name2, **metrics}
        )
        pbar.update(1)
    pbar.close()

    part_df = pd.DataFrame(records)
    with out_file.open("wb") as f:
        pickle.dump(part_df, f)
    return str(out_file)


@slurm_it(
    conda_env="multi-padlock-design",
    from_imports={"lib.nupack_heterodimer": "aggregate_nupack_selected_vs_all"},
)
def aggregate_nupack_selected_vs_all(
    output_dir: str,
    probe_df_path: str,
    sequence_col: str = "padlock",
    selected_indices: Optional[List[int]] = None,
    selected_gene_names: Optional[List[str]] = None,
    gene_col: str = "gene_name",
    dg_matrix_pickle: Optional[str] = None,
    percent_matrix_pickle: Optional[str] = None,
    combined_pickle: Optional[str] = None,
    overwrite: bool = True,
    use_seq_cache: bool = True,
    seq_cache_name: str = "probe_sequences_cache.pkl",
) -> str:
    """Aggregate per-probe selected-vs-all part files into row-subset vs all matrices.

    Rows correspond to the selected subset; columns correspond to all probes.
    """
    out_dir = Path(output_dir)
    if dg_matrix_pickle is None:
        dg_matrix_pickle = out_dir / "selected_vs_all_dg_matrix.pkl"
    else:
        dg_matrix_pickle = Path(dg_matrix_pickle)
    if percent_matrix_pickle is None:
        percent_matrix_pickle = out_dir / "selected_vs_all_percent_matrix.pkl"
    else:
        percent_matrix_pickle = Path(percent_matrix_pickle)
    if combined_pickle is None:
        combined_pickle = out_dir / "selected_vs_all_results.pkl"
    else:
        combined_pickle = Path(combined_pickle)

    if (
        dg_matrix_pickle.exists()
        and percent_matrix_pickle.exists()
        and combined_pickle.exists()
        and not overwrite
    ):
        return str(combined_pickle)

    probe_df_path = Path(probe_df_path)
    cache_file = out_dir / seq_cache_name
    names: Optional[List[str]] = None
    genes: Optional[List[str]] = None
    if use_seq_cache and cache_file.exists():
        try:
            with cache_file.open("rb") as f:
                cache_obj = pickle.load(f)
            if cache_obj.get("sequence_col") == sequence_col:
                names = cache_obj["names"]
        except Exception:
            names = None

    if names is None:
        if probe_df_path.suffix.lower() in {".pkl", ".pickle"}:
            df = pd.read_pickle(probe_df_path)
        else:
            df = pd.read_csv(probe_df_path)
        if sequence_col not in df.columns:
            raise ValueError(f"Sequence column '{sequence_col}' not found")
        name_col = None
        for c in ["name", "probe_name", "id", "gene", "padlock_name", "_autoname"]:
            if c in df.columns:
                name_col = c
                break
        if name_col is None:
            name_col = "_autoname"
            df[name_col] = [f"probe_{i}" for i in range(len(df))]
        names = df[name_col].astype(str).tolist()
        if (
            selected_indices is None
            and selected_gene_names is not None
            and gene_col in df.columns
        ):
            genes = df[gene_col].astype(str).tolist()

    N = len(names)

    # Determine selected indices if not given
    if selected_indices is None:
        if selected_gene_names is None:
            raise ValueError("Provide either selected_indices or selected_gene_names")
        if genes is None:
            if probe_df_path.suffix.lower() in {".pkl", ".pickle"}:
                df = pd.read_pickle(probe_df_path)
            else:
                df = pd.read_csv(probe_df_path)
            if gene_col not in df.columns:
                raise ValueError(f"Gene column '{gene_col}' not found")
            genes_series = df[gene_col].astype(str)
        else:
            genes_series = pd.Series(genes)
        gene_set = set(map(str, selected_gene_names))
        selected_indices = [i for i, g in enumerate(genes_series) if g in gene_set]

    selected_indices = list(sorted(set(selected_indices)))
    M = len(selected_indices)
    if M == 0:
        raise ValueError("No selected probes for aggregation")

    row_names = [names[i] for i in selected_indices]
    col_names = names

    dg_matrix = np.full((M, N), np.nan, dtype=float)
    percent_matrix = np.full((M, N), np.nan, dtype=float)

    # Map selected index to row position
    row_pos = {idx: r for r, idx in enumerate(selected_indices)}

    part_files = sorted(out_dir.glob("selected_vs_all_probe_*.pkl"))
    if not part_files:
        raise FileNotFoundError(
            "No per-probe part files found (selected_vs_all_probe_*.pkl)"
        )

    # Count rows for progress bar
    total_rows = 0
    for pf in part_files:
        try:
            with pf.open("rb") as f:
                part_df = pickle.load(f)
        except Exception:
            continue
        if isinstance(part_df, pd.DataFrame):
            total_rows += len(part_df)

    try:
        from tqdm.auto import tqdm
    except ImportError:

        def tqdm(x, **k):
            return x

    pbar = tqdm(
        total=total_rows,
        desc="aggregate selected-vs-all",
        unit="pair",
        ascii=True,
        disable=False,
        leave=False,
    )
    for pf in part_files:
        try:
            with pf.open("rb") as f:
                part_df = pickle.load(f)
        except Exception:
            continue
        if not isinstance(part_df, pd.DataFrame):
            continue
        # Each part file is a single i vs all j
        # Fill only when i belongs to selected set
        if "i" not in part_df.columns or "j" not in part_df.columns:
            continue
        i_values = part_df["i"].unique()
        if len(i_values) != 1:
            # Unexpected, but handle first value
            i_val = int(i_values[0])
        else:
            i_val = int(i_values[0])
        if i_val not in row_pos:
            pbar.update(len(part_df))
            continue
        r = row_pos[i_val]
        for row in part_df.itertuples():
            j = int(row.j)
            dg = getattr(row, "heterodimer_dG", math.nan)
            pct = getattr(row, "heterodimer_percent", math.nan)
            if 0 <= j < N:
                dg_matrix[r, j] = dg
                percent_matrix[r, j] = pct
            pbar.update(1)
    pbar.close()

    # Set diagonals (selected i where j==i) to 0.0 as convention
    for i in selected_indices:
        r = row_pos.get(i)
        if r is not None and 0 <= i < N:
            dg_matrix[r, i] = 0.0
            percent_matrix[r, i] = 0.0

    with open(dg_matrix_pickle, "wb") as f:
        pickle.dump(
            {
                "matrix": dg_matrix,
                "row_names": row_names,
                "col_names": col_names,
                "sequence_col": sequence_col,
            },
            f,
        )
    with open(percent_matrix_pickle, "wb") as f:
        pickle.dump(
            {
                "matrix": percent_matrix,
                "row_names": row_names,
                "col_names": col_names,
                "sequence_col": sequence_col,
            },
            f,
        )
    with open(combined_pickle, "wb") as f:
        pickle.dump(
            {
                "dg_matrix": dg_matrix,
                "percent_matrix": percent_matrix,
                "row_names": row_names,
                "col_names": col_names,
                "sequence_col": sequence_col,
            },
            f,
        )
    return str(combined_pickle)


def submit_selected_vs_all_jobs(
    probe_df_path: str,
    output_dir: str,
    slurm_folder: str,
    target_gene_names: List[str],
    gene_col: str = "gene_name",
    time: str = "12:00:00",
    mem: str = "8G",
    partition: Optional[str] = "ncpu",
    cpus_per_task: int = 1,
    sequence_col: str = "padlock",
    padlock_conc: float = 1e-9,
    dependency_aggregate_time: str = "01:00:00",
    dependency_aggregate_mem: str = "16G",
    dry_run: bool = False,
    use_seq_cache: bool = True,
    seq_cache_name: str = "probe_sequences_cache.pkl",
) -> tuple[list[str], Optional[str]]:
    """Submit one job per selected probe (selected-vs-all)
    plus dependent aggregator job.

    Returns the list of batch job IDs and the aggregator job ID.
    """
    probe_df_path_p = Path(probe_df_path)
    if probe_df_path_p.suffix.lower() in {".pkl", ".pickle"}:
        df = pd.read_pickle(probe_df_path_p)
    else:
        df = pd.read_csv(probe_df_path_p)
    if sequence_col not in df.columns:
        raise ValueError(f"Sequence column '{sequence_col}' not found")
    if gene_col not in df.columns:
        raise ValueError(f"Gene column '{gene_col}' not found")

    gene_set = set(map(str, target_gene_names))
    selected_indices = [
        int(i) for i in df.index[df[gene_col].astype(str).isin(gene_set)]
    ]
    if not selected_indices:
        raise ValueError(
            "No probes match the requested genes for selected-vs-all submission"
        )

    if dry_run:
        print(f"[DRY RUN] selected M={len(selected_indices)} of N={len(df)}")
        return [], None

    # Slurm options
    slurm_opts_batch = {"time": time, "mem": mem, "cpus-per-task": cpus_per_task}
    if partition:
        slurm_opts_batch["partition"] = partition
    slurm_opts_agg = {
        "time": dependency_aggregate_time,
        "mem": dependency_aggregate_mem,
        "cpus-per-task": 1,
    }
    if partition:
        slurm_opts_agg["partition"] = partition

    batch_job_ids: list[str] = []
    for idx in selected_indices:
        jid = run_nupack_selected_vs_all_for_probe(
            probe_df_path=str(probe_df_path_p),
            output_dir=str(output_dir),
            probe_index=idx,
            sequence_col=sequence_col,
            padlock_conc=padlock_conc,
            slurm_folder=slurm_folder,
            scripts_name=f"selected_vs_all_probe_{idx:06d}",
            slurm_options=slurm_opts_batch,
            use_slurm=True,
            use_seq_cache=use_seq_cache,
            seq_cache_name=seq_cache_name,
        )
        batch_job_ids.append(jid)

    aggregate_job_id = aggregate_nupack_selected_vs_all(
        output_dir=str(output_dir),
        probe_df_path=str(probe_df_path_p),
        sequence_col=sequence_col,
        selected_indices=selected_indices,
        selected_gene_names=list(gene_set),
        gene_col=gene_col,
        slurm_folder=slurm_folder,
        scripts_name="selected_vs_all_aggregate",
        job_dependency=batch_job_ids,
        slurm_options=slurm_opts_agg,
        use_slurm=True,
        use_seq_cache=use_seq_cache,
        seq_cache_name=seq_cache_name,
    )
    return batch_job_ids, aggregate_job_id


__all__ = [
    "submit_exhaustive_heterodimer_jobs",
    "run_nupack_heterodimer_batch",
    "aggregate_nupack_heterodimer_results",
    "run_nupack_selected_vs_all_for_probe",
    "aggregate_nupack_selected_vs_all",
    "submit_selected_vs_all_jobs",
]
