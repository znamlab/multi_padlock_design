"""Split a CSV file where the first column value becomes the filename.

Usage:
  python csv_splitter.py <input_csv> [--out OUTPUT_DIR]
Writes each row to <OUTPUT_DIR>/<first_col>.csv (single-row CSV).
If --out not provided, output directory is split_input/<input_stem>.
"""

import csv
import os
import sys
from pathlib import Path
from typing import Iterable

from multi_padlock_design.config import split_input


def write_row(out_dir: Path, row: Iterable[str]):
    cols = list(row)
    if not cols:
        return
    filename = out_dir / f"{cols[0]}.csv"
    with open(filename, "w", newline="") as newfile:
        csvwriter = csv.writer(newfile)
        csvwriter.writerow(cols)
    print(f"Created file: {filename}")


def split_csv(input_csv: str, output_dir: Path | None = None):
    if not os.path.isfile(input_csv):
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")
    input_path = Path(input_csv)
    if output_dir is None:
        output_dir = split_input / input_path.stem
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(input_csv, "r", newline="") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if row:
                write_row(output_dir, row)
    return output_dir


def parse_args(argv):
    if len(argv) < 2:
        print(__doc__)
        sys.exit(1)
    input_csv = None
    output_dir = None
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--out":
            if i + 1 >= len(argv):
                print("--out requires a path")
                sys.exit(1)
            output_dir = Path(argv[i + 1])
            i += 2
        elif arg.startswith("-"):
            print(f"Unknown option {arg}")
            sys.exit(1)
        else:
            if input_csv is None:
                input_csv = arg
            else:
                print("Unexpected extra positional argument")
                sys.exit(1)
            i += 1
    if input_csv is None:
        print("Input CSV required.")
        sys.exit(1)
    return input_csv, output_dir


def main():
    try:
        input_csv, output_dir = parse_args(sys.argv)
        out_dir = split_csv(input_csv, output_dir)
        print(f"All rows written under: {out_dir}")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
