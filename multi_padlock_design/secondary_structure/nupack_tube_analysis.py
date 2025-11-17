"""
NUPACK-based padlock analysis utilities extracted from the notebook.

Functions:
- analyze_padlock_df
- compute_capture_percentages
- analyze_isolated_padlocks
- analyze_padlock_bridging

Internal helpers are kept module-level to support multiprocessing pickling.
"""

from __future__ import annotations

import math
import re
from typing import Any, Dict, Iterable, Optional, Tuple

import pandas as pd
from nupack import Model, SetSpec, Strand, Tube, tube_analysis

from multi_padlock_design.secondary_structure import check_padlocks

# -----------------------------
# Naming helpers
# -----------------------------


def _extract_padlock_number(padlock_name: str) -> str:
    """Extract trailing digit group from padlock_name; return 'X' if none."""
    m = re.search(r"(\d+)(?!.*\d)", str(padlock_name))
    return m.group(1) if m else "X"


def _target_strand_name(gene_name: str, padlock_name: str) -> str:
    num = _extract_padlock_number(padlock_name)
    return f"{gene_name}_{num}"


# -----------------------------
# DataFrame-driven tube analysis
# -----------------------------


def analyze_padlock_df(
    df: pd.DataFrame,
    model: Optional[Model] = None,
    padlock_conc: float = 1e-9,
    target_conc: float = 1e-15,
    max_complex_size: int = 2,
    compute: Iterable[str] = ("pairs", "mfe"),
) -> tuple:
    """Build and analyze a single tube from a padlock/target dataframe.

    Target names are made unique per padlock as "gene_name_<digits>"; if no
    digits exist, 'X' is used.

    Args:
        df: DataFrame with columns {'padlock', 'padlock_name', 'gene_name',
            'target'}.
        model: NUPACK Model. If None, a default DNA model at 45 C (75 mM Na+,
            10 mM Mg2+) is created.
        padlock_conc: Molar concentration for padlock strands.
        target_conc: Molar concentration for reverse-complement target strands.
        max_complex_size: Maximum complex size allowed in SetSpec (>=2).
        compute: Iterable of NUPACK compute tasks (e.g., 'pairs', 'mfe', 'pfunc').

    Returns:
        Tuple[tube_result, tube, strands, df_with_names]:
            - tube_result: Result mapping from tube_analysis.
            - tube: The Tube that was analyzed.
            - strands: Dict[str, Strand] of all strands.
            - df_with_names: Copy of df with 'target_strand_name' added.

    Raises:
        ValueError: If required columns are missing from df.
    """
    required_cols = {"padlock", "padlock_name", "gene_name", "target"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Dataframe missing required columns: {missing}")

    if model is None:
        model = Model(material="dna", celsius=45, sodium=0.075, magnesium=0.01)

    df = df.copy()
    df["target_strand_name"] = [
        _target_strand_name(r["gene_name"], r["padlock_name"]) for _, r in df.iterrows()
    ]

    strands: Dict[str, Strand] = {}
    concs: Dict[Strand, float] = {}

    # Padlock strands
    for _, row in df.iterrows():
        seq = str(row["padlock"]).strip().upper()
        name = str(row["padlock_name"])
        if name in strands:
            suffix = 2
            new_name = f"{name}_{suffix}"
            while new_name in strands:
                suffix += 1
                new_name = f"{name}_{suffix}"
            name = new_name
        s = Strand(seq, name=name)
        strands[name] = s
        concs[s] = padlock_conc

    # Target (reverse complement) strands with unique names
    for _, row in df.iterrows():
        target_seq = str(row["target"]).strip().upper()
        name = row["target_strand_name"]
        rc_seq = check_padlocks.revcomp(target_seq)
        if name in strands:
            # Should not normally happen, but disambiguate
            suffix = 2
            new_name = f"{name}_{suffix}"
            while new_name in strands:
                suffix += 1
                new_name = f"{name}_{suffix}"
            name = new_name
            row["target_strand_name"] = name  # update copy
        s = Strand(rc_seq, name=name)
        strands[name] = s
        concs[s] = target_conc

    tube = Tube(
        strands=concs,
        complexes=SetSpec(max_size=max_complex_size),
        name="Padlock Tube",
    )
    tube_result = tube_analysis(tubes=[tube], compute=list(compute), model=model)
    return tube_result, tube, strands, df


def compute_capture_percentages(
    tube_result,
    tube: Tube,
    df: pd.DataFrame,
    fixed_target_conc: float = 1e-15,
) -> pd.DataFrame:
    """Compute capture percent for each padlock/target dimer in a tube.

    fraction = [padlock+target]_eq / fixed_target_conc; percent = fraction*100.

    Args:
        tube_result: Result mapping produced by tube_analysis for the tube.
        tube: The Tube used to extract complex concentrations.
        df: DataFrame with 'target_strand_name' (from analyze_padlock_df).
        fixed_target_conc: Normalization concentration for fraction/percent.

    Returns:
        DataFrame with columns:
            padlock_name, gene_name, target_strand_name, fixed_target_conc_M,
            padlock_target_complex_M, fraction, percent.

    Raises:
        ValueError: If 'target_strand_name' is missing from df.
    """
    if "target_strand_name" not in df.columns:
        raise ValueError(
            "DataFrame must contain 'target_strand_name' (run analyze_padlock_df first)"
        )

    tube_res = tube_result[tube]
    complex_concs = getattr(tube_res, "complex_concentrations", {})

    # Index dimer concentrations by ordered pair
    dimer_index: Dict[Tuple[str, str], float] = {}
    for cx, conc in complex_concs.items():
        try:
            names = tuple(s.name for s in cx.strands)
        except AttributeError:
            names = tuple(s.name for s in cx)
        if len(names) == 2:
            dimer_index[names] = conc

    def get_dimer(a: str, b: str) -> float:
        v = dimer_index.get((a, b))
        if v is not None:
            return v
        return dimer_index.get((b, a), math.nan)

    rows = []
    for _, r in df.iterrows():
        padlock = str(r["padlock_name"])
        target_name = str(r["target_strand_name"])
        gene = str(r["gene_name"])
        complex_m = get_dimer(padlock, target_name)
        if (not math.isnan(complex_m)) and fixed_target_conc > 0:
            fraction = complex_m / fixed_target_conc
            percent = fraction * 100
        else:
            fraction = math.nan
            percent = math.nan
        rows.append(
            {
                "padlock_name": padlock,
                "gene_name": gene,
                "target_strand_name": target_name,
                "fixed_target_conc_M": fixed_target_conc,
                "padlock_target_complex_M": complex_m,
                "fraction": fraction,
                "percent": percent,
            }
        )

    return pd.DataFrame(rows)


# -----------------------------
# Per-padlock isolated tube analysis (one tube per row)
# -----------------------------


def _analyze_isolated_worker(args: Tuple[Any, ...]) -> Dict[str, Any]:
    """Worker for a single padlock (separate process when parallel>1).

    Returns a dict row for the result DataFrame.
    """
    (
        idx,
        row,
        model_kwargs,
        padlock_conc,
        target_conc,
        max_complex_size,
        compute,
    ) = args

    # Import locally inside worker to be robust in subprocess context
    from nupack import Model, SetSpec, Strand, Tube, tube_analysis

    from ...lib import check_padlocks as _cp

    model = Model(**model_kwargs)
    padlock_seq = str(row["padlock"]).strip().upper()
    padlock_name = str(row["padlock_name"])
    gene_name = str(row["gene_name"])
    target_seq = str(row["target"]).strip().upper()

    try:
        target_name = _target_strand_name(gene_name, padlock_name)
    except Exception:
        target_name = f"{gene_name}_{idx}"

    padlock_strand = Strand(padlock_seq, name=padlock_name)
    target_rev = _cp.revcomp(target_seq)
    target_strand = Strand(target_rev, name=target_name)

    tube = Tube(
        strands={padlock_strand: padlock_conc, target_strand: target_conc},
        complexes=SetSpec(max_size=max_complex_size),
        name=f"Tube_{padlock_name}",
    )
    result_map = tube_analysis(tubes=[tube], compute=list(compute), model=model)

    # Get concentrations for complexes in this tube
    tube_res = result_map[tube]
    complex_concs = getattr(tube_res, "complex_concentrations", {})
    dimer_conc = math.nan
    for cx, conc in complex_concs.items():
        try:
            names = [s.name for s in cx.strands]
        except AttributeError:
            names = [s.name for s in cx]
        if len(names) == 2 and padlock_name in names and target_name in names:
            dimer_conc = conc
            break

    # Free energy: check common key string forms
    fe = math.nan
    expected1 = f"({padlock_name}+{target_name})"
    expected2 = f"({target_name}+{padlock_name})"
    for key in (expected1, expected2):
        if key in result_map:
            fe = getattr(result_map[key], "free_energy", math.nan)
            break

    if (not math.isnan(dimer_conc)) and target_conc > 0:
        fraction = dimer_conc / target_conc
        percent = fraction * 100
    else:
        fraction = math.nan
        percent = math.nan

    return {
        "padlock_name": padlock_name,
        "gene_name": gene_name,
        "target_strand_name": target_name,
        "fixed_target_conc_M": target_conc,
        "padlock_target_complex_M": dimer_conc,
        "fraction": fraction,
        "percent": percent,
        "free_energy_complex_kcal_per_mol": fe,
    }


def analyze_isolated_padlocks(
    df: pd.DataFrame,
    model: Optional[Model] = None,
    padlock_conc: float = 1e-9,
    target_conc: float = 1e-15,
    max_complex_size: int = 2,
    compute: Iterable[str] = ("pfunc", "mfe"),
    progress: bool = False,
    parallel: int = 1,
) -> pd.DataFrame:
    """Analyze each padlock in isolation with its target in separate tubes.

    Runs one tube per row (padlock + rc(target)). Determines 'pfunc'
    which gives free energy value; supports multiprocessing but is buggy.

    Args:
        df: DataFrame with columns {'padlock','padlock_name','gene_name','target'}.
        model: NUPACK Model. If None, uses a default DNA model (45 C, 75 mM Na+,
            10 mM Mg2+).
        padlock_conc: Molar concentration for the padlock in each tube.
        target_conc: Molar concentration for the target in each tube.
        max_complex_size: Maximum complex size (>=2).
        compute: Iterable of compute tasks (e.g., 'pfunc', 'mfe').
        progress: If True, show a progress bar.
        parallel: Number of processes to use (>1 enables multiprocessing).

    Returns:
        DataFrame with columns:
            padlock_name, gene_name, target_strand_name, fixed_target_conc_M,
            padlock_target_complex_M, fraction, percent,
            free_energy_complex_kcal_per_mol.

    Raises:
        ValueError: If required columns are missing from df.
    """
    if progress:
        try:
            from tqdm.auto import tqdm
        except ImportError:

            def tqdm(x, **k):
                return x

    else:

        def tqdm(x, **k):
            return x

    required_cols = {"padlock", "padlock_name", "gene_name", "target"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    if model is None:
        model = Model(material="dna", celsius=45, sodium=0.075, magnesium=0.01)

    if "pfunc" not in compute:
        compute = tuple(["pfunc", *compute])

    # Sequential path
    model_kwargs = dict(material="dna", celsius=45, sodium=0.075, magnesium=0.01)
    if parallel <= 1:
        iter_rows = df.iterrows()
        iter_rows = tqdm(iter_rows, total=len(df), desc="Analyzing padlocks")
        records = []
        for idx, row in iter_rows:
            rec = _analyze_isolated_worker(
                (
                    idx,
                    row,
                    model_kwargs,
                    padlock_conc,
                    target_conc,
                    max_complex_size,
                    compute,
                )
            )
            records.append(rec)
        return pd.DataFrame(records)

    # Parallel path
    from multiprocessing import Pool

    jobs = [
        (idx, row, model_kwargs, padlock_conc, target_conc, max_complex_size, compute)
        for idx, row in df.iterrows()
    ]

    records = []
    with Pool(processes=parallel) as pool:
        if progress:
            from tqdm.auto import tqdm as _tqdm

            for rec in _tqdm(
                pool.imap_unordered(_analyze_isolated_worker, jobs),
                total=len(jobs),
                desc=f"Analyzing padlocks (parallel x{parallel})",
            ):
                records.append(rec)
        else:
            for rec in pool.imap_unordered(_analyze_isolated_worker, jobs):
                records.append(rec)

    # Preserve original order
    return (
        pd.DataFrame(records)
        .set_index("padlock_name")
        .reindex([str(r["padlock_name"]) for _, r in df.iterrows()])
        .reset_index(drop=False)
    )


# -----------------------------
# Bridging analysis (both padlock arms bound simultaneously)
# -----------------------------


def _split_revcomp_target(full_target_seq: str) -> Tuple[str, str]:
    """Return two halves of the reverse complement of the provided target sequence.

    The split is length//2 for the first half and the remainder for the second, so that:
        len(half1) + len(half2) == len(rc_full)
    This assumes the padlock arms bind adjacent halves of the target.
    Adjust here if you have precise arm lengths.
    """
    rc = check_padlocks.revcomp(full_target_seq.upper().strip())
    mid = len(rc) // 2
    return rc[:mid], rc[mid:]


def _bridging_worker(args: Tuple[Any, ...]) -> Dict[str, Any]:
    """Worker performing tube_analysis for a single padlock bridging scenario."""
    (
        idx,
        row,
        model_kwargs,
        padlock_conc,
        target_conc,
        max_complex_size,
        compute,
    ) = args

    model = Model(**model_kwargs)
    padlock_seq = str(row["padlock"]).strip().upper()
    padlock_name = str(row["padlock_name"])
    gene_name = str(row["gene_name"])
    target_seq = str(row["target"]).strip().upper()

    try:
        base_target_name = _target_strand_name(gene_name, padlock_name)
    except Exception:
        base_target_name = f"{gene_name}_{idx}"

    half1_seq, half2_seq = _split_revcomp_target(target_seq)
    half1_name = f"{base_target_name}_H1"
    half2_name = f"{base_target_name}_H2"

    padlock_strand = Strand(padlock_seq, name=padlock_name)
    half1_strand = Strand(half1_seq, name=half1_name)
    half2_strand = Strand(half2_seq, name=half2_name)

    # Need complexes up to size 3 to allow padlock bridging both halves simultaneously
    if max_complex_size < 3:
        max_complex_size = 3

    tube = Tube(
        strands={
            padlock_strand: padlock_conc,
            half1_strand: target_conc,
            half2_strand: target_conc,
        },
        complexes=SetSpec(max_size=max_complex_size),
        name=f"Bridge_{padlock_name}",
    )

    result_map = tube_analysis(tubes=[tube], compute=list(compute), model=model)
    tube_res = result_map[tube]
    complex_concs = getattr(tube_res, "complex_concentrations", {})

    dimer1 = math.nan  # padlock + half1
    dimer2 = math.nan  # padlock + half2
    triple = math.nan  # padlock + half1 + half2 (bridged)

    # Iterate complexes to find relevant species
    for cx, conc in complex_concs.items():
        try:
            names = [s.name for s in cx.strands]
        except AttributeError:
            names = [s.name for s in cx]
        if padlock_name in names:
            if len(names) == 2:
                if half1_name in names:
                    dimer1 = conc
                elif half2_name in names:
                    dimer2 = conc
            elif len(names) == 3 and half1_name in names and half2_name in names:
                triple = conc

    # Free energy for triple complex if present (any ordering)
    triple_fe = math.nan
    if not math.isnan(triple):
        expected_forms = [
            f"({padlock_name}+{half1_name}+{half2_name})",
            f"({padlock_name}+{half2_name}+{half1_name})",
            f"({half1_name}+{padlock_name}+{half2_name})",
            f"({half1_name}+{half2_name}+{padlock_name})",
            f"({half2_name}+{padlock_name}+{half1_name})",
            f"({half2_name}+{half1_name}+{padlock_name})",
        ]
        for key in expected_forms:
            if key in result_map:
                triple_fe = getattr(result_map[key], "free_energy", math.nan)
                break

    # Fractions: triple complex relative to initial target concentration
    # (single full target)
    # Assumption: Each half introduced at target_conc but we still normalize
    # to target_conc
    if target_conc > 0 and not math.isnan(triple):
        triple_fraction = triple / target_conc
        triple_percent = triple_fraction * 100
    else:
        triple_fraction = math.nan
        triple_percent = math.nan

    # Optional: fractions of individual dimers (any binding) relative to target_conc
    frac_half1_any = (
        (dimer1 + (triple if not math.isnan(triple) else 0)) / target_conc
        if target_conc > 0 and (not math.isnan(dimer1) or not math.isnan(triple))
        else math.nan
    )
    frac_half2_any = (
        (dimer2 + (triple if not math.isnan(triple) else 0)) / target_conc
        if target_conc > 0 and (not math.isnan(dimer2) or not math.isnan(triple))
        else math.nan
    )

    return {
        "padlock_name": padlock_name,
        "gene_name": gene_name,
        "target_half1_name": half1_name,
        "target_half2_name": half2_name,
        "half1_length": len(half1_seq),
        "half2_length": len(half2_seq),
        "padlock_conc_M": padlock_conc,
        "half_initial_conc_M": target_conc,
        "triple_complex_M": triple,
        "triple_fraction": triple_fraction,
        "triple_percent": triple_percent,
        "half1_dimer_M": dimer1,
        "half2_dimer_M": dimer2,
        "half1_any_fraction": frac_half1_any,
        "half2_any_fraction": frac_half2_any,
        "free_energy_triple_kcal_per_mol": triple_fe,
    }


def analyze_padlock_bridging(
    df: pd.DataFrame,
    model: Optional[Model] = None,
    padlock_conc: float = 1e-9,
    target_conc: float = 1e-15,
    max_complex_size: int = 3,
    compute: Iterable[str] = ("pfunc",),
    progress: bool = False,
    parallel: int = 1,
) -> pd.DataFrame:
    """Analyze simultaneous binding of both padlock arms (bridging).

    For each row, rc(target) is split into halves (H1, H2); a tube with
    {padlock, H1, H2} is analyzed allowing size >= 3 complexes. The concentration
    of the triple complex (padlock+H1+H2) is used to compute bridging fractions.

    Args:
        df: DataFrame with columns {'padlock','padlock_name','gene_name','target'}.
        model: NUPACK Model. If None, uses a default DNA model (45 C, 75 mM Na+,
            10 mM Mg2+).
        padlock_conc: Molar concentration for the padlock strand.
        target_conc: Molar concentration assigned to each half target strand.
        max_complex_size: Maximum complex size (forced to >=3 for bridging).
        compute: Iterable of compute tasks (ensures 'pfunc').
        progress: If True, show a progress bar.
        parallel: Number of processes to use (>1 enables multiprocessing).

    Returns:
        DataFrame with columns:
            padlock_name, gene_name, target_half1_name, target_half2_name,
            half1_length, half2_length, padlock_conc_M, half_initial_conc_M,
            triple_complex_M, triple_fraction, triple_percent,
            half1_dimer_M, half2_dimer_M, half1_any_fraction, half2_any_fraction,
            free_energy_triple_kcal_per_mol.

    Raises:
        ValueError: If required columns are missing from df.
    """
    if progress:
        try:
            from tqdm.auto import tqdm
        except ImportError:

            def tqdm(x, **k):
                return x

    else:

        def tqdm(x, **k):
            return x

    required = {"padlock", "padlock_name", "gene_name", "target"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if model is None:
        model = Model(material="dna", celsius=45, sodium=0.075, magnesium=0.01)

    if "pfunc" not in compute:
        compute = tuple(["pfunc", *compute])

    if max_complex_size < 3:
        max_complex_size = 3

    model_kwargs = dict(material="dna", celsius=45, sodium=0.075, magnesium=0.01)

    if parallel <= 1:
        iterator = df.iterrows()
        iterator = tqdm(iterator, total=len(df), desc="Bridging analysis")
        records = []
        for idx, row in iterator:
            rec = _bridging_worker(
                (
                    idx,
                    row,
                    model_kwargs,
                    padlock_conc,
                    target_conc,
                    max_complex_size,
                    compute,
                )
            )
            records.append(rec)
        return pd.DataFrame(records)

    # Parallel path
    from multiprocessing import Pool

    jobs = [
        (idx, row, model_kwargs, padlock_conc, target_conc, max_complex_size, compute)
        for idx, row in df.iterrows()
    ]
    records = []
    with Pool(processes=parallel) as pool:
        if progress:
            from tqdm.auto import tqdm as _tqdm

            for rec in _tqdm(
                pool.imap_unordered(_bridging_worker, jobs),
                total=len(jobs),
                desc=f"Bridging (parallel x{parallel})",
            ):
                records.append(rec)
        else:
            for rec in pool.imap_unordered(_bridging_worker, jobs):
                records.append(rec)

    # Restore original order
    ordered = (
        pd.DataFrame(records)
        .set_index("padlock_name")
        .reindex([str(r["padlock_name"]) for _, r in df.iterrows()])
        .reset_index(drop=False)
    )
    return ordered


__all__ = [
    "analyze_padlock_df",
    "compute_capture_percentages",
    "analyze_isolated_padlocks",
    "analyze_padlock_bridging",
]
