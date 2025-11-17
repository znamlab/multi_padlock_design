"""
Hybridisation probe design pipeline utilities.

Functions here generate random probe sequences, submit BLAST
searches via parallel_blast_query.sh, parse BLAST outputs, compute
nearest-neighbour Tm against subject hits (cseq), and summarize keep/reject.

Intended to be called from the hyb_probe_design.ipynb notebook.
"""

from __future__ import annotations

import random
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

from multi_padlock_design.blast import readblast

# ----------------------------- RNG & basic utils -----------------------------
_rng = random.Random()


def set_seed(seed: Optional[int]) -> None:
	"""Seed the internal RNG used for sequence generation."""
	global _rng
	_rng = random.Random(int(seed)) if seed is not None else random.Random()


DNA = "ACGT"
GC_SET = set("GC")


def gc_content(seq: str) -> float:
	seq = seq.upper()
	if not seq:
		return 0.0
	gc = sum(1 for b in seq if b in GC_SET)
	return gc / len(seq)


def has_gc_clamps(seq: str) -> bool:
	seq = seq.upper()
	return (
		len(seq) >= 4
		and set(seq[:2]).issubset(GC_SET)
		and set(seq[-2:]).issubset(GC_SET)
	)


def rand_seq_with_clamps(length: int) -> str:
	"""Generate a random DNA sequence with GC clamps at both ends."""
	prefix = "".join(_rng.choice("GC") for _ in range(2))
	suffix = "".join(_rng.choice("GC") for _ in range(2))
	n_mid = max(0, length - 4)
	middle = "".join(_rng.choice(DNA) for _ in range(n_mid))
	return prefix + middle + suffix


# ---------------------------- Sequence generation ---------------------------
def generate_candidates(
	n: int,
	min_len: int,
	max_len: int,
	gc_min: float,
	gc_max: float,
	require_gc_clamps: bool,
	tm_min: float,
	tm_max: float,
	max_attempts: int = 200000,
) -> pd.DataFrame:
	"""Generate up to n sequences satisfying GC%, length, GC-clamps and Tm constraints.

	Returns a DataFrame with: padlock_name, target, length, gc_frac, tm_self, gc_clamps
	"""
	rows = []
	attempts = 0
	while len(rows) < n and attempts < max_attempts:
		attempts += 1
		L = _rng.randint(min_len, max_len)
		if require_gc_clamps:
			seq = rand_seq_with_clamps(L)
			if not has_gc_clamps(seq):
				continue
		else:
			seq = "".join(_rng.choice(DNA) for _ in range(L))
		gc = gc_content(seq)
		if not (gc_min <= gc <= gc_max):
			continue
		try:
			tm_self = float(readblast.calc_tm_NN(seq=seq))
		except Exception:
			continue
		if not (tm_min <= tm_self <= tm_max):
			continue
		idx = len(rows) + 1
		rows.append(
			{
				"padlock_name": f"probe_{idx}",
				"target": seq,
				"length": L,
				"gc_frac": gc,
				"tm_self": tm_self,
				"gc_clamps": True if require_gc_clamps else has_gc_clamps(seq),
			}
		)
	if len(rows) < n:
		print(
			"Warning: generated only {} sequences after {} attempts; "
			"consider relaxing constraints".format(len(rows), attempts)
		)
	return pd.DataFrame(rows)


def write_probe_csv_for_blast(
	df: pd.DataFrame, work_dir: Path, gene_tag: str
) -> Path:
	"""Write candidates to a CSV consumed by padlocks_to_fasta.py.

	Columns: gene_name, padlock_name, target
	"""
	out_csv = Path(work_dir) / f"{gene_tag}.csv"
	tmp = df.copy()
	tmp.insert(0, "gene_name", gene_tag)
	tmp = tmp[["gene_name", "padlock_name", "target"]]
	out_csv.parent.mkdir(parents=True, exist_ok=True)
	tmp.to_csv(out_csv, index=False)
	return out_csv


# --------------------------------- BLAST I/O --------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PAR_SCRIPT = REPO_ROOT / "padlock_checking" / "parallel_blast_query.sh"


def submit_parallel_blast(
	csv_path: Path,
	work_dir: Path,
	par_script: Optional[Path] = None,
	sbatch: bool = True,
) -> Optional[str]:
	"""Submit the parallel BLAST job(s) via the project script using sbatch.

	Returns sbatch output string if submitted, else None.
	"""
	par = Path(par_script) if par_script else DEFAULT_PAR_SCRIPT
	if not par.is_file():
		raise FileNotFoundError(f"parallel_blast_query.sh not found: {par}")
	if not sbatch:
		print("sbatch disabled; run manually:")
		print(" sbatch", par, csv_path, work_dir)
		return None
	cmd = ["sbatch", str(par), str(csv_path), str(Path(work_dir))]
	print("Submitting:", " ".join(cmd))
	out = subprocess.check_output(cmd, text=True).strip()
	print(out)
	return out


def expected_blast_outputs(work_dir: Path) -> List[Path]:
	"""Return expected *_blast.out files (using fasta_list.txt when present)."""
	work_dir = Path(work_dir)
	fasta_list = work_dir / "fasta_list.txt"
	outs: List[Path] = []
	if fasta_list.exists():
		with open(fasta_list) as fh:
			for line in fh:
				p = Path(line.strip())
				if p.name.endswith(".fasta"):
					outs.append(p.with_name(p.stem + "_blast.out"))
	else:
		outs = sorted(work_dir.glob("*_blast.out"))
	return outs


def wait_for_blast_outputs(
	work_dir: Path,
	poll_seconds: int = 20,
	timeout_seconds: Optional[int] = None,
) -> List[Path]:
	"""Poll for BLAST outputs to appear and be non-empty.

	Returns the expected outputs when detected (non-empty files).
	"""
	outs = expected_blast_outputs(work_dir)
	if not outs:
		print("No expected outputs detected; waiting for any *_blast.out in work_dir…")
	waited = 0
	start = time.time()
	while True:
		if outs and all(p.exists() and p.stat().st_size > 0 for p in outs):
			print("BLAST jobs appear complete.")
			return outs
		# fallback: if none listed, finish when any exist
		if not outs:
			any_out = list(Path(work_dir).glob("*_blast.out"))
			if any_out:
				print("Detected some BLAST outputs; proceeding.")
				return any_out
		if timeout_seconds is not None and (time.time() - start) > timeout_seconds:
			print("Timeout reached while waiting for BLAST outputs.")
			return outs
		time.sleep(poll_seconds)
		waited += poll_seconds
		print(f"...still waiting ({waited}s)")


# ----------------------------- Parsing & scoring ----------------------------
def parse_blast_outputs(work_dir: Path) -> pd.DataFrame:
	"""Parse all *_blast.out in work_dir into one DataFrame.
	"""
	from glob import glob

	work = Path(work_dir)
	blast_files = sorted(glob(str(work / "*_blast.out")))
	if not blast_files:
		print("No BLAST outputs found in", work)
		return pd.DataFrame()

	dfs = []
	for bf in blast_files:
		try:
			dfb = pd.read_csv(bf, header=None)
		except Exception:
			continue
		if dfb.empty:
			continue
		# outfmt "10 std qseq sseq" -> first 14 cols
		if dfb.shape[1] < 14:
			continue
		dfb = dfb.iloc[:, :14]
		dfb.columns = [
			"length",
			"mismatch",
			"gapopen",
			"qstart",
			"qend",
			"sstart",
			"send",
			"evalue",
			"bitscore",
			"qseq",
			"sseq",
		]
		dfb["blast_file"] = bf
		dfs.append(dfb)
	return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def compute_tm_hit_cseq(parsed: pd.DataFrame) -> pd.DataFrame:
	"""Add tm_hit_cseq column (Tm of qseq vs sseq using readblast.calc_tm_NN)."""
	if parsed is None or parsed.empty:
		return parsed

	def _tm(row) -> float:
		try:
			return float(readblast.calc_tm_NN(seq=row["qseq"], cseq=row["sseq"]))
		except Exception:
			return np.nan

	out = parsed.copy()
	out["tm_hit_cseq"] = out.apply(_tm, axis=1)
	return out


def summarize_candidates(
	df_candidates: pd.DataFrame,
	parsed_with_tm: pd.DataFrame,
	blast_tm_cutoff: float,
	work_dir: Path,
	gene_tag: str,
) -> Tuple[pd.DataFrame, Path]:
	"""Merge candidates with max per-hit Tm and decide keep/reject.

	Also writes CSV summary to work_dir. Returns (df_out, csv_path).
	"""
	if parsed_with_tm is not None and not parsed_with_tm.empty:
		parsed_with_tm = parsed_with_tm.copy()
		parsed_with_tm["padlock_name"] = parsed_with_tm["qaccver"].astype(str)
		max_tm_per_query = (
			parsed_with_tm.groupby("padlock_name")["tm_hit_cseq"].max().rename("max_blast_tm_cseq")
		)
	else:
		max_tm_per_query = pd.Series(dtype=float)

	df_out = df_candidates.copy()
	df_out["padlock_name"] = df_out["padlock_name"].astype(str)
	df_out = df_out.merge(
		max_tm_per_query, how="left", left_on="padlock_name", right_index=True
	)
	df_out["max_blast_tm_cseq"] = df_out["max_blast_tm_cseq"].fillna(-np.inf)
	df_out["kept"] = df_out["max_blast_tm_cseq"] <= float(blast_tm_cutoff)
	df_out["reason"] = np.where(
		df_out["kept"], "no hit above cutoff", "rejected: high-Tm hit"
	)

	out_table = Path(work_dir) / f"{gene_tag}_probe_screen.csv"
	df_out.to_csv(out_table, index=False)
	return df_out, out_table


# ----------------------------- Convenience runner ---------------------------
def make_job_tag(prefix: str = "rand_probes") -> str:
	return f"{prefix}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
