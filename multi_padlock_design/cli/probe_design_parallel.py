from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Sequence

from multi_padlock_design import config
from multi_padlock_design.io import csv_splitter, fasta_splitter


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    description = (
        "Split a FASTA or CSV file into per-record inputs and submit a SLURM job for "
        "each chunk using scripts/probe_design.sh."
    )
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_file",
        help="Path to the .fa/.fasta/.csv file to split before submission.",
    )
    parser.add_argument(
        "--log-root",
        type=Path,
        help=(
            "Override where SLURM stdout files are written. "
            "Defaults to <repo>/logs/slurm_logs/<input_stem>."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print sbatch commands without submitting jobs.",
    )
    return parser.parse_args(argv)


def split_inputs(input_path: Path, input_type: str) -> Path:
    destination = config.split_input / input_path.stem
    destination.mkdir(parents=True, exist_ok=True)

    if input_type == "fasta":
        fasta_splitter.split_fasta(str(input_path), output_dir=destination)
    else:
        csv_splitter.split_csv(str(input_path), output_dir=destination)

    return destination


def collect_split_files(out_dir: Path, suffixes: Iterable[str]) -> List[Path]:
    lowered = {suffix.lower() for suffix in suffixes}
    return sorted(
        p for p in out_dir.iterdir() if p.is_file() and p.suffix.lower() in lowered
    )


def build_export_env(
    parent_base: str,
    input_type: str,
    chunk_path: Path,
) -> str:
    defaults = config.probe_design_defaults
    env_pairs = {
        "INPUT": str(chunk_path.resolve()),
        "PARENT": parent_base,
        "INPUT_TYPE": input_type,
        "SPECIES": defaults["species"],
        "ARM_LEN": defaults["arm_len"],
        "TOTAL_LEN": defaults["total_len"],
        "INTERVAL": defaults["interval"],
        "TM_LOW": defaults["tm_low"],
        "TM_HIGH": defaults["tm_high"],
        "N_PROBES": defaults["n_probes"],
    }
    return "--export=" + ",".join(f"{k}={v}" for k, v in env_pairs.items())


def submit_jobs(
    repo_root: Path,
    input_type: str,
    chunk_files: Sequence[Path],
    log_dir: Path,
    dry_run: bool = False,
) -> None:
    script_path = repo_root / "scripts" / "probe_design.sh"
    if not script_path.is_file():
        raise FileNotFoundError(
            f"Missing probe_design.sh under {repo_root / 'scripts'}"
        )

    parent_base = chunk_files[0].parent.name
    for chunk in chunk_files:
        base_name = chunk.stem
        export_arg = build_export_env(parent_base, input_type, chunk)
        log_file = log_dir / f"{parent_base}_{base_name}.out"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            "sbatch",
            export_arg,
            f"--output={log_file}",
            str(script_path),
        ]
        if dry_run:
            printable = " ".join(shlex.quote(part) for part in cmd)
            print(f"DRY RUN: {printable}")
            continue
        subprocess.run(cmd, check=True)
        print(f"Submitted job for {chunk}")


def determine_input_type(path: Path) -> str:
    ext = path.suffix.lower().lstrip(".")
    if ext in {"fa", "fasta"}:
        return "fasta"
    if ext == "csv":
        return "csv"
    raise ValueError("Input file must have .fa, .fasta, or .csv extension.")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    input_path = Path(args.input_file).expanduser().resolve()
    if not input_path.is_file():
        print(f"Input file not found: {input_path}", file=sys.stderr)
        return 1

    try:
        input_type = determine_input_type(input_path)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    repo_root = config.REPO_ROOT

    split_dir = split_inputs(input_path, input_type)
    suffixes = [".fasta"] if input_type == "fasta" else [".csv"]
    chunk_files = collect_split_files(split_dir, suffixes)

    if not chunk_files:
        print(f"No split files found under {split_dir}")
        return 0

    if args.log_root is not None:
        log_root = args.log_root.expanduser().resolve()
    else:
        log_root = repo_root / "logs" / "slurm_logs"
    log_dir = log_root / split_dir.name
    log_dir.mkdir(parents=True, exist_ok=True)

    try:
        submit_jobs(repo_root, input_type, chunk_files, log_dir, dry_run=args.dry_run)
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"Error submitting jobs: {exc}", file=sys.stderr)
        return 1

    print(f"Queued {len(chunk_files)} jobs for input '{input_path.name}'.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
