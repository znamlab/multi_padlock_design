import random
import re

import pandas as pd

import multi_padlock_design.blast.readblast as readblast


def revcomp(seq):
    complement = str.maketrans("ATCG", "TAGC")
    return seq.translate(complement)[::-1]


def read_candidates(outdir, panel_name, probefile="3.AllSpecificTargets_*.csv"):
    """
    Read and concatenate all "3.AllSpecificTargets_*.csv" files from subdirectories of
    the given panel directory.

    Args:
        outdir (Path): Base output directory containing the panel subdirectory.
        panel_name (str): Name of the panel subdirectory.

    Returns:
        pd.DataFrame: Concatenated dataframe with columns
        ['target', 'Tm', 'startpos', 'gene'].
    """
    all_dfs = []
    panel_dir = outdir / panel_name

    # Loop over subdirectories
    for subdir in panel_dir.iterdir():
        if subdir.is_dir():
            print(f"Processing {subdir.name}")

            # Collect "3.AllSpecificTargets_*.csv"
            target_files = list(subdir.glob(probefile))
            assert (
                len(target_files) == 1
            ), f"Expected exactly one target file in {subdir},"
            f" found {len(target_files)}"
            tf = target_files[0]

            # Load and tag with gene name (= subdir name)
            # remove FASTA-like header lines first
            with open(tf) as f:
                clean_lines = [line for line in f if not line.startswith(">")]
                clean_lines = [smart_split(line.strip()) for line in clean_lines]

            header, *data = clean_lines
            df = pd.DataFrame(data, columns=header)
            df["gene"] = subdir.name
            if len(df) < 1:
                print(f"Warning: No targets found in {subdir.name}, skipping.")
                continue
            all_dfs.append(df)

        # Concatenate all into one big dataframe
        big_df = pd.concat(all_dfs, ignore_index=True)

    return big_df


def smart_split(line):
    # split on commas not inside [ ... ]
    return re.split(r",(?![^[]*\])", line)


def generate_random_dna_sequences(
    min_barcode, length, num_sequences, for_overlap, genes=None
):
    """
    Build the primer-padlock overlapping sequences for SNAIL probes
    """

    bases = ["A", "T", "C", "G"]
    sequences = []
    reverse_complements = []
    Tms = []

    while len(sequences) < num_sequences:
        sequence = "".join(random.choice(bases) for _ in range(length))

        # Check condition for overlap oligos
        if for_overlap and sequence[3:6] in [
            "TAA",
            "TCA",
            "TGA",
            "TTA",
            "ATT",
            "AGT",
            "ACT",
            "AAT",
        ]:  # ligation won't happen when reverse complement
            # is TXA In 6:9 position for overlap
            continue

        # No homopolymers
        if (
            ("AAA" in sequence)
            or ("TTT" in sequence)
            or ("CCC" in sequence)
            or ("GGG" in sequence)
        ):
            continue

        # check Tm?
        Tm = readblast.calc_tm_NN(sequence)
        if Tm < 25 or Tm > 35:
            continue

        if all(
            count_diff(seq[:min_barcode], sequence[:min_barcode]) >= 2
            for seq in sequences
        ) and all(count_diff(seq, sequence) >= 4 for seq in sequences):
            sequences.append(sequence)
            complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
            reverse_seq = sequence[::-1]
            reverse_complement_seq = "".join(complement[base] for base in reverse_seq)
            reverse_complements.append(reverse_complement_seq)
            Tms.append(Tm)

    # Build dataframe
    df = pd.DataFrame(
        {"Sequence": sequences, "Reverse Complement": reverse_complements, "Tms": Tms}
    ).drop_duplicates()

    # Handle gene names
    if genes is None:
        genes = ["no_gene"] * len(df)
    else:
        # pad or truncate to match df length
        if len(genes) < len(df):
            genes = genes + ["no_gene"] * (len(df) - len(genes))
        elif len(genes) > len(df):
            genes = genes[: len(df)]

    df["Overlap Primer Sequence"] = genes
    df["ID"] = range(len(df))  # unique IDs

    return df


def count_diff(seq1, seq2):
    return sum(base1 != base2 for base1, base2 in zip(seq1, seq2))
