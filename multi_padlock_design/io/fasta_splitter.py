import os
import sys
from pathlib import Path
from typing import List

from Bio import SeqIO

import multi_padlock_design.config as config

"""
Split a multi-record FASTA file into per-gene FASTA files.
  - Gene name taken as the first whitespace-delimited token of the record.description
  - Sequence uppercased
  - Header simplified to just the gene name (id, name, description all set)
  - Records for the same gene are appended to that gene's FASTA file
Also produces a `gene_names.txt` listing unique gene names (in encounter order)
unless disabled.

Usage:
    python fasta_splitter.py <input_fasta> [--out OUTPUT_DIR] [--no-list]
Examples:
  python fasta_splitter.py olfr_consensus_cds_3utr.fa
  python fasta_splitter.py olfr_consensus_cds_3utr.fa /tmp/output_fastas
  python fasta_splitter.py olfr_consensus_cds_3utr.fa --no-list
  python fasta_splitter.py olfr_consensus_cds_3utr.fa /tmp/output_fastas --no-list
"""


def extract_gene_name(description: str) -> str:
    """Return the first token from a FASTA description line."""
    return description.split()[0] if description else "UNKNOWN"


def split_fasta(
    input_fasta: str,
    output_dir: Path | None = None,
    write_list: bool = True,
) -> List[str]:
    if not os.path.isfile(input_fasta):
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

    input_path = Path(input_fasta)
    if output_dir is None:
        # Derive directory name from input file stem under split_input root
        output_dir = config.split_input / input_path.stem
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_gene_names: List[str] = []

    # Iterate through records and write out
    for record in SeqIO.parse(input_fasta, "fasta"):
        gene_name = extract_gene_name(record.description)
        raw_gene_names.append(gene_name)

        # Uppercase sequence
        record.seq = record.seq.upper()

        # Simplify header
        record.id = gene_name
        record.name = gene_name
        record.description = gene_name

        out_path = output_dir / f"{gene_name}.fasta"
        # Append mode to allow multiple records per gene
        with open(out_path, "a") as handle:
            SeqIO.write(record, handle, "fasta")

    # Deduplicate preserving order
    seen = set()
    gene_names = [g for g in raw_gene_names if not (g in seen or seen.add(g))]

    if write_list:
        list_path = output_dir / "gene_names.txt"
        with open(list_path, "w") as fh:
            fh.write("\n".join(gene_names))
    else:
        list_path = ""

    print(
        f"Processed {len(raw_gene_names)} records; "
        f"{len(gene_names)} unique gene names."
    )
    print(f"Per-gene FASTA files written to: {output_dir}")
    if write_list:
        print(f"Gene list written to: {list_path}")

    return gene_names


def parse_args(argv: List[str]):
    """Parse CLI arguments."""
    if len(argv) < 2:
        print(__doc__)
        sys.exit(1)

    input_fasta = None
    output_dir: Path | None = None
    write_list = True

    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == "--no-list":
            write_list = False
            i += 1
        elif arg == "--out":
            if i + 1 >= len(argv):
                print("--out requires a path")
                sys.exit(1)
            output_dir = Path(argv[i + 1])
            i += 2
        elif arg.startswith("-"):
            print(f"Unknown option: {arg}")
            sys.exit(1)
        else:
            if input_fasta is None:
                input_fasta = arg
            else:
                print("Unexpected extra positional argument.")
                sys.exit(1)
            i += 1

    if input_fasta is None:
        print("Input FASTA required.")
        sys.exit(1)

    if output_dir is None:
        # Default: directory containing the input FASTA
        print("Output directory is required")

    return input_fasta, output_dir, write_list


def main():
    try:
        input_fasta, output_dir, write_list = parse_args(sys.argv)
        split_fasta(input_fasta, output_dir, write_list=write_list)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
