import argparse
from pathlib import Path

import pandas as pd


def write_fasta_from_padlocks(padlock_file_path: str, output_dir: str):
    """
    Read padlocks CSV (expects columns: gene_name, padlock_name, target)
    and write one FASTA per gene plus a fasta_list.txt.
    """
    padlock_file = Path(padlock_file_path)
    if not padlock_file.is_file():
        raise FileNotFoundError(f"Padlock file not found: {padlock_file}")

    padlocks = pd.read_csv(padlock_file, header=0)
    # Tolerate unnamed index column
    padlocks = padlocks.loc[:, ~padlocks.columns.str.contains("^Unnamed")]
    required = {"gene_name", "padlock_name", "target"}
    missing = required - set(padlocks.columns)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for gene in padlocks.gene_name.dropna().unique():
        tmp = padlocks.loc[
            padlocks["gene_name"] == gene, ["target", "gene_name", "padlock_name"]
        ]
        fasta_path = out_dir / f"{gene}_query.fasta"
        print(f"Writing {fasta_path}", flush=True)
        with open(fasta_path, "w") as fh:
            for padlock in tmp.padlock_name:
                seq = tmp.loc[tmp["padlock_name"] == padlock, "target"].iloc[0]
                fh.write(f">{padlock}\n{seq}\n")

    fasta_paths = sorted(out_dir.glob("*.fasta"))
    list_path = out_dir / "fasta_list.txt"
    with open(list_path, "w") as fh:
        for p in fasta_paths:
            fh.write(str(p) + "\n")
    print(f"Wrote FASTA list: {list_path}")
    return list_path


def parse_args():
    ap = argparse.ArgumentParser(
        description="Generate per-gene FASTA files from padlock CSV."
    )
    ap.add_argument("--input", "-i", required=True, help="Padlocks CSV file")
    ap.add_argument(
        "--output", "-o", required=True, help="Output directory for FASTA files"
    )
    return ap.parse_args()


def main():
    args = parse_args()
    write_fasta_from_padlocks(args.input, args.output)


if __name__ == "__main__":
    main()
