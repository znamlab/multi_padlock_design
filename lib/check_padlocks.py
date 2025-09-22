import os
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
from Bio.Seq import Seq
import subprocess
from pathlib import Path

from lib.readblast import has_gap_or_mismatch, split_arms, fill_gaps, calc_tm_NN
from lib import formatrefseq
import config


def find_off_targets(gene, blast_query_path, armlength=20):
    print(f"Finding off-targets for {gene}", flush=True)
    file = blast_query_path / f"{gene}_query_blast.out"

    df = pd.read_csv(file, header=None)

    # Exclude lines with predicted transcripts
    df = df.loc[~(df[1].str.contains("XR_") | df[1].str.contains("XM_"))]

    # Create a new column 'homology_candidate_hit' based on the given conditions

    # Exclude hits with less than 50% coverage and less than 80% homology

    df["homology_candidate_hit"] = (
        # over 80% homology
        (df[2] > 80)
        &
        # over 50% coverage
        (df[3] > 2 * armlength * 0.5)
        &
        # hit covers +/- 5bp from ligation site
        (df[6] < armlength - 4)
        & (df[7] > armlength + 5)
    )
    failed_match = False
    if getattr(config, "reference_transcriptome", "").lower() == "refseq":
        variants = find_variants(gene)
    elif getattr(config, "reference_transcriptome", "").lower() == "ensembl":
        if getattr(config, "annotation_file", None):
            # Convert Olfr gene name to gene symbol
            try:
                anno_df = pd.read_csv(config.annotation_file)
                if {"gene_name", "gene_symbol"}.issubset(anno_df.columns):
                    match = anno_df.loc[anno_df["gene_name"] == gene, "gene_symbol"]
                    if (
                        not match.empty
                        and pd.notna(match.values[0])
                        and str(match.values[0]).strip()
                    ):
                        gene = match.values[0]
                    else:
                        warnings.warn(
                            f"Gene '{gene}' not found in annotation_file; using original gene name."
                        )
                        failed_match = True
                else:
                    warnings.warn(
                        "annotation_file missing required columns 'gene_name' and 'gene_symbol'; using original gene name."
                    )
            except Exception as e:
                warnings.warn(
                    f"Failed to map gene using annotation_file ({e}); using original gene name."
                )
        variants = find_variants(gene)
    # Exclude known variants
    # df = df.loc[~(df[1].isin(variants))]

    df["gene"] = gene

    # Write off-targets to file
    df.to_csv(blast_query_path / f"{gene}_off_targets.out")

    if failed_match:
        return gene
    else:
        return None


def find_variants(gene):
    """Find transcript/accession variants corresponding to a gene.
    For Olfr genes:
      - Rows containing '|' are split (take the part before the first '|').
      - For those rows, if an annotation file is provided (expects columns: gene_name,gene_symbol),
        attempt to replace the split name with its gene_symbol.
      - Other rows remain unchanged.
    """
    ref_path = os.path.join(config.BASE_DIR, config.reference_transcriptome)
    acronyms = pd.read_csv(
        os.path.join(ref_path, f"{config.species}.acronymheaders.txt"), header=None
    )
    acronyms.columns = ["gene_name"]

    # Identify rows containing '|'
    mask_pipe = acronyms["gene_name"].str.contains("|", regex=False)
    pipe_indices = acronyms.index[mask_pipe]
    if not pipe_indices.empty:
        base_names = acronyms.loc[pipe_indices, "gene_name"].str.split("|").str[0]

        # Prepare replacements: default to base name
        replacements = dict(zip(pipe_indices, base_names))

        # Optional gene_name -> gene_symbol mapping
        if getattr(config, "annotation_file", None):
            try:
                anno_df = pd.read_csv(config.annotation_file)
                if {"gene_name", "gene_symbol"}.issubset(anno_df.columns):
                    mapping = dict(zip(anno_df["gene_name"], anno_df["gene_symbol"]))
                    for idx in pipe_indices:
                        base = replacements[idx]
                        gene_symbol = mapping.get(base)
                        if gene_symbol:  # replace only if a symbol exists
                            replacements[idx] = gene_symbol
                # else silently fall back to base names
            except Exception:
                # On any failure keep base names
                pass

        # Apply replacements
        for idx, new_name in replacements.items():
            acronyms.at[idx, "gene_name"] = new_name

    # Find indices matching the (possibly already symbol) gene
    indices = acronyms.index[acronyms.gene_name == gene].to_list()

    accession_numbers = []
    with open(
        os.path.join(ref_path, f"{config.species}.selectedheaders.txt"), "r"
    ) as f:
        lines = f.readlines()
        for i in indices:
            accession_numbers.append(lines[int(i)].strip())

    accession_numbers = format_variants(accession_numbers)
    return accession_numbers


def format_variants(accession_numbers):
    # Remove the > character and take everything before the first whitespace character
    accession_numbers = [accession.replace(">", "") for accession in accession_numbers]
    accession_numbers = [accession.split(" ")[0] for accession in accession_numbers]
    return accession_numbers


# Define the process_row function at the top level
def process_row(row, armlength=20, tm_threshold=37, precomputed_variants=None):
    query_seq = row.qseq
    subject_seq = row.sseq
    ligation_site = armlength + 1
    query_start = row.qstart
    hit = row.subject  # Get the subject (hit)

    # Use the precomputed variants for this query
    variants = precomputed_variants[row.gene]

    valid_probe = False
    query_left = None
    query_right = None
    tm_left = None
    tm_right = None
    tm_no_mismatch_left = None
    tm_no_mismatch_right = None
    specific = None  # Default specific value

    # Split the sequences into left and right arms using query-derived split for both
    (
        query_left,
        query_right,
        subject_left,
        subject_right,
    ) = split_arms(query_seq, subject_seq, ligation_site, query_start)

    # Check for gaps or mismatches near the ligation site
    if not has_gap_or_mismatch(query_seq, subject_seq, ligation_site, query_start):
        ligation_site_missmatch = False

        # Fill gaps and calculate Tm for left arm
        if query_left:
            query_left, subject_left = fill_gaps(query_left, subject_left)
            tm_left = calc_tm_NN(seq=query_left, cseq=subject_left)
            tm_no_mismatch_left = calc_tm_NN(seq=query_left)

        # Fill gaps and calculate Tm for right arm
        if query_right:
            query_right, subject_right = fill_gaps(query_right, subject_right)
            tm_right = calc_tm_NN(seq=query_right, cseq=subject_right)
            tm_no_mismatch_right = calc_tm_NN(seq=query_right)

        # Check if both arms have Tm values greater than or equal to the threshold
        if (tm_left is not None and tm_left >= tm_threshold) and (
            tm_right is not None and tm_right >= tm_threshold
        ):
            valid_probe = True

        else:
            valid_probe = False  # Tm is below the threshold, invalid probe
    else:
        ligation_site_missmatch = (
            True  # There is a gap or mismatch near the ligation site
        )
    # Check if variants are provided
    if len(variants):
        # If the hit is not in the variants, this is a non-specific hit
        if hit not in variants:
            specific = False
        else:
            # Otherwise, it's specific
            specific = True
    else:
        if row.blast_target == row.gene:
            specific = True
        else:
            specific = False
    return {
        "query_left": query_left,
        "query_right": query_right,
        "subject_left": subject_left,
        "subject_right": subject_right,
        "tm_left_NN": tm_left,
        "tm_right_NN": tm_right,
        "tm_no_mismatch_left_NN": tm_no_mismatch_left,
        "tm_no_mismatch_right_NN": tm_no_mismatch_right,
        "valid_probe_NN": valid_probe,
        "ligation_site_missmatch": ligation_site_missmatch,
        "specific": specific,  # Add specific column to the result
    }


# --- New global cache for fast bulk variant lookup ---
_GENE_VARIANT_INDEX = None


def _build_gene_variant_index():
    """
    Build and return a dict: gene_name (or mapped gene_symbol) -> list of accession variants.
    Replicates the transformation logic from find_variants but does it once for all genes.
    """
    ref_path = os.path.join(config.BASE_DIR, config.reference_transcriptome)
    acronym_path = os.path.join(ref_path, f"{config.species}.acronymheaders.txt")
    selected_headers_path = os.path.join(
        ref_path, f"{config.species}.selectedheaders.txt"
    )
    print("building")
    # Load acronym headers
    acronyms_df = pd.read_csv(acronym_path, header=None, names=["gene_name"])

    # Process '|' lines (same logic as find_variants)
    mask_pipe = acronyms_df["gene_name"].str.contains("|", regex=False)
    if mask_pipe.any():
        base_names = acronyms_df.loc[mask_pipe, "gene_name"].str.split("|").str[0]
        acronyms_df.loc[mask_pipe, "gene_name"] = base_names

        # Optional annotation mapping (gene_name -> gene_symbol)
        if getattr(config, "annotation_file", None):
            try:
                anno_df = pd.read_csv(config.annotation_file)
                if {"gene_name", "gene_symbol"}.issubset(anno_df.columns):
                    mapping = dict(zip(anno_df["gene_name"], anno_df["gene_symbol"]))
                    acronyms_df.loc[mask_pipe, "gene_name"] = acronyms_df.loc[
                        mask_pipe, "gene_name"
                    ].map(lambda x: mapping.get(x, x))
            except Exception:
                pass  # Fail silently, keep base names

    # Load selected headers once
    with open(selected_headers_path, "r") as f:
        selected_headers = [line.rstrip("\n") for line in f]

    # Build index
    variant_index = {}
    for idx, gene in enumerate(acronyms_df["gene_name"]):
        header = selected_headers[idx]
        # format like format_variants would (strip '>' and take token before first space)
        acc = header.replace(">", "").split(" ")[0]
        variant_index.setdefault(gene, []).append(acc)

    return variant_index


# Precompute the variants for each unique query
def precompute_variants(df):
    """
    Optimized version:
      - Builds a global gene->variants index once.
      - Performs O(1) lookups per unique query (no repeated file I/O).
    """
    global _GENE_VARIANT_INDEX
    if _GENE_VARIANT_INDEX is None:
        _GENE_VARIANT_INDEX = _build_gene_variant_index()
    print("test")
    unique_queries = df["gene"].unique()
    precomputed_variants = {}
    for query in tqdm(unique_queries, desc="Precomputing variants"):
        # Fallback to original find_variants if gene not in bulk index (edge cases)
        variants = _GENE_VARIANT_INDEX.get(query)
        if variants is None:
            variants = find_variants(query)
        precomputed_variants[query] = variants
    return precomputed_variants


# Vectorized processing without multiprocessing, using precomputed variants
def process_dataframe(df, armlength, tm_threshold, precomputed_variants):
    # Precompute variants for each unique query

    results = []

    # tqdm for showing progress
    for row in tqdm(df.itertuples(), total=df.shape[0], desc="Processing rows"):
        result = process_row(row, armlength, tm_threshold, precomputed_variants)
        results.append(result)

    # Create a new DataFrame from the results list
    results_df = pd.DataFrame(results)

    df = df.reset_index(drop=True)
    results_df = results_df.reset_index(drop=True)

    # Concatenate the results with the original DataFrame
    df = pd.concat([df, results_df], axis=1)

    return df


# Cached data
cached_data = {}


def loaddb(species, config):
    """Load formatted RefSeq header and sequence files"""
    if species == "mouse":
        s = 0
    elif species == "human":
        s = 1
    else:
        raise ValueError(f"Unknown species: {species}")

    if species not in cached_data:
        try:
            # Determine whether to use RefSeq or Ensembl
            if config.reference_transcriptome == "refseq":
                fastadir = os.path.join(config.BASE_DIR, config.reference_transcriptome)
                fasta_filenum = config.fasta_filenum_refseq
                fasta_pre_suffix = config.fasta_pre_suffix_refseq
                # create trascriptome database if not existing
                if not os.path.isfile(fastadir + "/" + species + ".acronymheaders.txt"):
                    print("Processing fasta database files..")
                    formatrefseq.fastadb(
                        fastadir, fasta_filenum, fasta_pre_suffix, species
                    )
            elif config.reference_transcriptome == "ensembl":
                fastadir = os.path.join(config.BASE_DIR, config.reference_transcriptome)
                files = [config.cdna_file]
                if hasattr(config, "extra_files") and config.extra_files:
                    files += config.extra_files
                if not os.path.isfile(fastadir + "/" + species + ".acronymheaders.txt"):
                    print("Processing fasta database files..")
                    formatrefseq.fastadb(fastadir, files, species)

            # load database files
            with open(fastadir + "/" + species + ".acronymheaders.txt", "r") as f:
                Acronyms = [line.rstrip("\n") for line in f]
            with open(fastadir + "/" + species + ".selectedheaders.txt", "r") as f:
                Headers = [line.rstrip("\n") for line in f]
            with open(fastadir + "/" + species + ".selectedseqs.txt", "r") as f:
                Seq = [line.rstrip("\n") for line in f]
        except FileNotFoundError:
            print("Cannot load fasta database due to mismatch in species.")

        cached_data[species] = (Acronyms, Headers, Seq)
    return cached_data[species]


def find_genes_from_variants(variants, species, config):
    """
    Reverse lookup: given a list of variant identifiers (accession IDs, raw
    acronym/header entries, or pipe-delimited Ensembl-style lines like
      >Olfr1307|OTTMUSG00000015117|OTTMUST00000035855
    return their corresponding gene symbol(s) (or processed gene names) using
    the same transformation rules applied in find_variants:
      - For acronym lines containing '|', only the part before the first '|'
        is kept.
      - If annotation_file (with columns gene_name,gene_symbol) is provided,
        base names are mapped to gene_symbol when possible.
    The function tolerates leading '>' and ignores text after the first
    whitespace in headers. Output: dict {original_input: [gene_symbol] or None}.
    """
    Acronyms, Headers, _ = loaddb(species, config)

    # Build DataFrame for Acronyms to apply same transformations.

    acr_df = pd.DataFrame({"raw": Acronyms})
    # Extract base gene name (split at first '|' if present)
    has_pipe = acr_df["raw"].str.contains("|", regex=False)
    acr_df["base"] = acr_df["raw"]
    acr_df.loc[has_pipe, "base"] = acr_df.loc[has_pipe, "raw"].str.split("|").str[0]

    # Optional annotation mapping
    mapping = {}
    if getattr(config, "annotation_file", None):
        try:
            anno_df = pd.read_csv(config.annotation_file)
            if {"gene_name", "gene_symbol"}.issubset(anno_df.columns):
                mapping = dict(zip(anno_df["gene_name"], anno_df["gene_symbol"]))
        except Exception:
            mapping = {}

    # Apply mapping to base -> symbol where available
    acr_df["processed"] = acr_df["base"].map(lambda x: mapping.get(x, x))

    # Build accession mapping aligned by index with Headers (same ordering as in loaddb)
    # Accession token: strip '>' then take text before first whitespace
    def _header_to_accession(h):
        h = h.lstrip(">")
        return h.split(" ")[0]

    accessions = [_header_to_accession(h) for h in Headers]

    # Maps:
    accession_to_gene = {}
    raw_acronym_to_gene = {}
    for idx, proc in enumerate(acr_df["processed"]):
        accession_to_gene[accessions[idx]] = proc
        raw_acronym = acr_df.at[idx, "raw"]
        raw_acronym_to_gene[raw_acronym] = proc
        # Also allow processed form itself as key (self-mapping)
        raw_acronym_to_gene.setdefault(proc, proc)

    result = {}
    for original in variants:
        v = str(original).strip()
        if not v:
            result[original] = None
            continue

        # Normalize: remove leading '>' and isolate first whitespace token
        v_norm = v.lstrip(">").split()[0]

        gene_candidates = None

        # 1. Direct accession match
        if v_norm in accession_to_gene:
            gene_candidates = [accession_to_gene[v_norm]]
        # 2. Raw acronym match
        elif v_norm in raw_acronym_to_gene:
            gene_candidates = [raw_acronym_to_gene[v_norm]]
        else:
            # 3. Pipe-delimited header provided as variant: take base before first '|'
            if "|" in v_norm:
                base = v_norm.split("|", 1)[0]
            else:
                base = v_norm
            # Map base through annotation mapping if possible
            processed = mapping.get(base, base)
            if processed in raw_acronym_to_gene:
                gene_candidates = [raw_acronym_to_gene[processed]]
            else:
                # Accept processed itself if appears anywhere in processed list
                if processed in acr_df["processed"].values:
                    gene_candidates = [processed]

        result[original] = gene_candidates

    return result


dna_IMM1 = {
    "AG/TT": (1.0, 0.9),
    "AT/TG": (-2.5, -8.3),
    "CG/GT": (-4.1, -11.7),
    "CT/GG": (-2.8, -8.0),
    "GG/CT": (3.3, 10.4),
    "GG/TT": (5.8, 16.3),
    "GT/CG": (-4.4, -12.3),
    "GT/TG": (4.1, 9.5),
    "TG/AT": (-0.1, -1.7),
    "TG/GT": (-1.4, -6.2),
    "TT/AG": (-1.3, -5.3),
    "AA/TG": (-0.6, -2.3),
    "AG/TA": (-0.7, -2.3),
    "CA/GG": (-0.7, -2.3),
    "CG/GA": (-4.0, -13.2),
    "GA/CG": (-0.6, -1.0),
    "GG/CA": (0.5, 3.2),
    "TA/AG": (0.7, 0.7),
    "TG/AA": (3.0, 7.4),
    "AC/TT": (0.7, 0.2),
    "AT/TC": (-1.2, -6.2),
    "CC/GT": (-0.8, -4.5),
    "CT/GC": (-1.5, -6.1),
    "GC/CT": (2.3, 5.4),
    "GT/CC": (5.2, 13.5),
    "TC/AT": (1.2, 0.7),
    "TT/AC": (1.0, 0.7),
    "AA/TC": (2.3, 4.6),
    "AC/TA": (5.3, 14.6),
    "CA/GC": (1.9, 3.7),
    "CC/GA": (0.6, -0.6),
    "GA/CC": (5.2, 14.2),
    "GC/CA": (-0.7, -3.8),
    "TA/AC": (3.4, 8.0),
    "TC/AA": (7.6, 20.2),
    "AA/TA": (1.2, 1.7),
    "CA/GA": (-0.9, -4.2),
    "GA/CA": (-2.9, -9.8),
    "TA/AA": (4.7, 12.9),
    "AC/TC": (0.0, -4.4),
    "CC/GC": (-1.5, -7.2),
    "GC/CC": (3.6, 8.9),
    "TC/AC": (6.1, 16.4),
    "AG/TG": (-3.1, -9.5),
    "CG/GG": (-4.9, -15.3),
    "GG/CG": (-6.0, -15.8),
    "TG/AG": (1.6, 3.6),
    "AT/TT": (-2.7, -10.8),
    "CT/GT": (-5.0, -15.8),
    "GT/CT": (-2.2, -8.4),
    "TT/AT": (0.2, -1.5),
    "AI/TC": (-8.9, -25.5),
    "TI/AC": (-5.9, -17.4),
    "AC/TI": (-8.8, -25.4),
    "TC/AI": (-4.9, -13.9),
    "CI/GC": (-5.4, -13.7),
    "GI/CC": (-6.8, -19.1),
    "CC/GI": (-8.3, -23.8),
    "GC/CI": (-5.0, -12.6),
    "AI/TA": (-8.3, -25.0),
    "TI/AA": (-3.4, -11.2),
    "AA/TI": (-0.7, -2.6),
    "TA/AI": (-1.3, -4.6),
    "CI/GA": (2.6, 8.9),
    "GI/CA": (-7.8, -21.1),
    "CA/GI": (-7.0, -20.0),
    "GA/CI": (-7.6, -20.2),
    "AI/TT": (0.49, -0.7),
    "TI/AT": (-6.5, -22.0),
    "AT/TI": (-5.6, -18.7),
    "TT/AI": (-0.8, -4.3),
    "CI/GT": (-1.0, -2.4),
    "GI/CT": (-3.5, -10.6),
    "CT/GI": (0.1, -1.0),
    "GT/CI": (-4.3, -12.1),
    "AI/TG": (-4.9, -15.8),
    "TI/AG": (-1.9, -8.5),
    "AG/TI": (0.1, -1.8),
    "TG/AI": (1.0, 1.0),
    "CI/GG": (7.1, 21.3),
    "GI/CG": (-1.1, -3.2),
    "CG/GI": (5.8, 16.9),
    "GG/CI": (-7.6, -22.0),
    "AI/TI": (-3.3, -11.9),
    "TI/AI": (0.1, -2.3),
    "CI/GI": (1.3, 3.0),
    "GI/CI": (-0.5, -1.3),
    "AA/AA": (3.8, 6.3),
    "AA/AC": (3.9, 6.3),
    "AA/AG": (3.8, 6.3),
    "AA/CA": (3.9, 6.6),
    "AA/CC": (4.1, 6.6),
    "AA/CG": (4.4, 6.6),
    "AA/GA": (3.3, 4.7),
    "AA/GC": (3.4, 4.7),
    "AA/GG": (3.3, 4.7),
    "AC/AA": (4.3, 7.5),
    "AC/AC": (4.4, 7.5),
    "AC/AT": (4.9, 7.5),
    "AC/CA": (4.5, 7.8),
    "AC/CC": (4.7, 7.8),
    "AC/CT": (4.6, 7.8),
    "AC/GA": (3.8, 5.9),
    "AC/GC": (3.9, 5.9),
    "AC/GT": (4.3, 5.9),
    "AG/AA": (4.1, 7.0),
    "AG/AG": (4.1, 7.0),
    "AG/AT": (4.7, 7.0),
    "AG/CA": (4.3, 7.3),
    "AG/CG": (4.8, 7.3),
    "AG/CT": (4.4, 7.3),
    "AG/GA": (3.6, 5.4),
    "AG/GG": (3.6, 5.4),
    "AG/GT": (4.1, 5.4),
    "AT/AC": (4.6, 6.6),
    "AT/AG": (4.5, 6.6),
    "AT/AT": (5.1, 6.6),
    "AT/CC": (4.9, 6.9),
    "AT/CG": (5.2, 6.9),
    "AT/CT": (4.8, 6.9),
    "AT/GC": (4.2, 5.0),
    "AT/GG": (4.1, 5.0),
    "AT/GT": (4.6, 5.0),
    "CA/AA": (4.0, 6.3),
    "CA/AC": (4.1, 6.3),
    "CA/AG": (4.0, 6.3),
    "CA/CA": (4.1, 6.6),
    "CA/CC": (4.3, 6.6),
    "CA/CG": (4.6, 6.6),
    "CA/TA": (2.0, -1.8),
    "CA/TC": (1.6, -1.8),
    "CA/TG": (2.0, -1.8),
    "CC/AA": (4.5, 7.5),
    "CC/AC": (4.6, 7.5),
    "CC/AT": (5.1, 7.5),
    "CC/CA": (4.7, 7.8),
    "CC/CC": (4.9, 7.8),
    "CC/CT": (4.8, 7.8),
    "CC/TA": (2.6, -0.6),
    "CC/TC": (2.2, -0.6),
    "CC/TT": (2.0, -0.6),
    "CG/AA": (4.3, 7.0),
    "CG/AG": (4.3, 7.0),
    "CG/AT": (4.9, 7.0),
    "CG/CA": (4.5, 7.3),
    "CG/CG": (5.0, 7.3),
    "CG/CT": (4.6, 7.3),
    "CG/TA": (2.4, -1.1),
    "CG/TG": (2.4, -1.1),
    "CG/TT": (1.8, -1.1),
    "CT/AC": (4.2, 6.6),
    "CT/AG": (4.1, 6.6),
    "CT/AT": (4.7, 6.6),
    "CT/CC": (4.5, 6.9),
    "CT/CG": (4.8, 6.9),
    "CT/CT": (4.4, 6.9),
    "CT/TC": (1.8, -1.5),
    "CT/TG": (2.2, -1.5),
    "CT/TT": (1.6, -1.5),
    "GA/AA": (3.9, 6.3),
    "GA/AC": (4.0, 6.3),
    "GA/AG": (3.9, 6.3),
    "GA/GA": (3.4, 4.7),
    "GA/GC": (3.5, 4.7),
    "GA/GG": (3.4, 4.7),
    "GA/TA": (1.9, -1.8),
    "GA/TC": (1.5, -1.8),
    "GA/TG": (1.9, -1.8),
    "GC/AA": (4.6, 7.5),
    "GC/AC": (4.7, 7.5),
    "GC/AT": (5.2, 7.5),
    "GC/GA": (4.1, 5.9),
    "GC/GC": (4.2, 5.9),
    "GC/GT": (4.6, 5.9),
    "GC/TA": (2.7, -0.6),
    "GC/TC": (2.3, -0.6),
    "GC/TT": (2.1, -0.6),
    "GG/AA": (4.1, 7.0),
    "GG/AG": (4.1, 7.0),
    "GG/AT": (4.7, 7.0),
    "GG/GA": (3.6, 5.4),
    "GG/GG": (3.6, 5.4),
    "GG/GT": (4.1, 5.4),
    "GG/TA": (2.2, -1.1),
    "GG/TG": (2.2, -1.1),
    "GG/TT": (1.6, -1.1),
    "GT/AC": (4.6, 6.6),
    "GT/AG": (4.5, 6.6),
    "GT/AT": (5.1, 6.6),
    "GT/GC": (4.2, 5.0),
    "GT/GG": (4.1, 5.0),
    "GT/GT": (4.6, 5.0),
    "GT/TC": (2.2, -1.5),
    "GT/TG": (2.6, -1.5),
    "GT/TT": (2.0, -1.5),
    "TA/CA": (4.6, 6.6),
    "TA/CC": (4.8, 6.6),
    "TA/CG": (5.1, 6.6),
    "TA/GA": (4.0, 4.7),
    "TA/GC": (4.1, 4.7),
    "TA/GG": (4.0, 4.7),
    "TA/TA": (2.5, -1.8),
    "TA/TC": (2.1, -1.8),
    "TA/TG": (2.5, -1.8),
    "TC/CA": (4.6, 7.8),
    "TC/CC": (4.8, 7.8),
    "TC/CT": (4.7, 7.8),
    "TC/GA": (3.9, 5.9),
    "TC/GC": (4.0, 5.9),
    "TC/GT": (4.4, 5.9),
    "TC/TA": (2.5, -0.6),
    "TC/TC": (2.1, -0.6),
    "TC/TT": (1.9, -0.6),
    "TG/CA": (4.9, 7.3),
    "TG/CG": (5.4, 7.3),
    "TG/CT": (5.0, 7.3),
    "TG/GA": (4.2, 5.4),
    "TG/GG": (4.2, 5.4),
    "TG/GT": (4.7, 5.4),
    "TG/TA": (2.8, -1.1),
    "TG/TG": (2.8, -1.1),
    "TG/TT": (2.2, -1.1),
    "TT/CC": (4.3, 6.9),
    "TT/CG": (4.6, 6.9),
    "TT/CT": (4.2, 6.9),
    "TT/GC": (3.6, 5.0),
    "TT/GG": (3.5, 5.0),
    "TT/GT": (4.0, 5.0),
    "TT/TC": (1.6, -1.5),
    "TT/TG": (2.0, -1.5),
    "TT/TT": (1.4, -1.5),
}

DNA_TMM1 = {
    "AA/TA": (-3.1, -7.8),
    "TA/AA": (-2.5, -6.3),
    "CA/GA": (-4.3, -10.7),
    "GA/CA": (-8.0, -22.5),
    "AC/TC": (-0.1, 0.5),
    "TC/AC": (-0.7, -1.3),
    "CC/GC": (-2.1, -5.1),
    "GC/CC": (-3.9, -10.6),
    "AG/TG": (-1.1, -2.1),
    "TG/AG": (-1.1, -2.7),
    "CG/GG": (-3.8, -9.5),
    "GG/CG": (-0.7, -19.2),
    "AT/TT": (-2.4, -6.5),
    "TT/AT": (-3.2, -8.9),
    "CT/GT": (-6.1, -16.9),
    "GT/CT": (-7.4, -21.2),
    "AA/TC": (-1.6, -4.0),
    "AC/TA": (-1.8, -3.8),
    "CA/GC": (-2.6, -5.9),
    "CC/GA": (-2.7, -6.0),
    "GA/CC": (-5.0, -13.8),
    "GC/CA": (-3.2, -7.1),
    "TA/AC": (-2.3, -5.9),
    "TC/AA": (-2.7, -7.0),
    "AC/TT": (-0.9, -1.7),
    "AT/TC": (-2.3, -6.3),
    "CC/GT": (-3.2, -8.0),
    "CT/GC": (-3.9, -10.6),
    "GC/CT": (-4.9, -13.5),
    "GT/CC": (-3.0, -7.8),
    "TC/AT": (-2.5, -6.3),
    "TT/AC": (-0.7, -1.2),
    "AA/TG": (-1.9, -4.4),
    "AG/TA": (-2.5, -5.9),
    "CA/GG": (-3.9, -9.6),
    "CG/GA": (-6.0, -15.5),
    "GA/CG": (-4.3, -11.1),
    "GG/CA": (-4.6, -11.4),
    "TA/AG": (-2.0, -4.7),
    "TG/AA": (-2.4, -5.8),
    "AG/TT": (-3.2, -8.7),
    "AT/TG": (-3.5, -9.4),
    "CG/GT": (-3.8, -9.0),
    "CT/GG": (-6.6, -18.7),
    "GG/CT": (-5.7, -15.9),
    "GT/CG": (-5.9, -16.1),
    "TG/AT": (-3.9, -10.5),
    "TT/AG": (-3.6, -9.8),
    "AA/AA": (3.8, 6.3),
    "AA/AC": (3.9, 6.3),
    "AA/AG": (3.8, 6.3),
    "AA/CA": (3.9, 6.6),
    "AA/CC": (4.1, 6.6),
    "AA/CG": (4.4, 6.6),
    "AA/GA": (3.3, 4.7),
    "AA/GC": (3.4, 4.7),
    "AA/GG": (3.3, 4.7),
    "AC/AA": (4.3, 7.5),
    "AC/AC": (4.4, 7.5),
    "AC/AT": (4.9, 7.5),
    "AC/CA": (4.5, 7.8),
    "AC/CC": (4.7, 7.8),
    "AC/CT": (4.6, 7.8),
    "AC/GA": (3.8, 5.9),
    "AC/GC": (3.9, 5.9),
    "AC/GT": (4.3, 5.9),
    "AG/AA": (4.1, 7.0),
    "AG/AG": (4.1, 7.0),
    "AG/AT": (4.7, 7.0),
    "AG/CA": (4.3, 7.3),
    "AG/CG": (4.8, 7.3),
    "AG/CT": (4.4, 7.3),
    "AG/GA": (3.6, 5.4),
    "AG/GG": (3.6, 5.4),
    "AG/GT": (4.1, 5.4),
    "AT/AC": (4.6, 6.6),
    "AT/AG": (4.5, 6.6),
    "AT/AT": (5.1, 6.6),
    "AT/CC": (4.9, 6.9),
    "AT/CG": (5.2, 6.9),
    "AT/CT": (4.8, 6.9),
    "AT/GC": (4.2, 5.0),
    "AT/GG": (4.1, 5.0),
    "AT/GT": (4.6, 5.0),
    "CA/AA": (4.0, 6.3),
    "CA/AC": (4.1, 6.3),
    "CA/AG": (4.0, 6.3),
    "CA/CA": (4.1, 6.6),
    "CA/CC": (4.3, 6.6),
    "CA/CG": (4.6, 6.6),
    "CA/TA": (2.0, -1.8),
    "CA/TC": (1.6, -1.8),
    "CA/TG": (2.0, -1.8),
    "CC/AA": (4.5, 7.5),
    "CC/AC": (4.6, 7.5),
    "CC/AT": (5.1, 7.5),
    "CC/CA": (4.7, 7.8),
    "CC/CC": (4.9, 7.8),
    "CC/CT": (4.8, 7.8),
    "CC/TA": (2.6, -0.6),
    "CC/TC": (2.2, -0.6),
    "CC/TT": (2.0, -0.6),
    "CG/AA": (4.3, 7.0),
    "CG/AG": (4.3, 7.0),
    "CG/AT": (4.9, 7.0),
    "CG/CA": (4.5, 7.3),
    "CG/CG": (5.0, 7.3),
    "CG/CT": (4.6, 7.3),
    "CG/TA": (2.4, -1.1),
    "CG/TG": (2.4, -1.1),
    "CG/TT": (1.8, -1.1),
    "CT/AC": (4.2, 6.6),
    "CT/AG": (4.1, 6.6),
    "CT/AT": (4.7, 6.6),
    "CT/CC": (4.5, 6.9),
    "CT/CG": (4.8, 6.9),
    "CT/CT": (4.4, 6.9),
    "CT/TC": (1.8, -1.5),
    "CT/TG": (2.2, -1.5),
    "CT/TT": (1.6, -1.5),
    "GA/AA": (3.9, 6.3),
    "GA/AC": (4.0, 6.3),
    "GA/AG": (3.9, 6.3),
    "GA/GA": (3.4, 4.7),
    "GA/GC": (3.5, 4.7),
    "GA/GG": (3.4, 4.7),
    "GA/TA": (1.9, -1.8),
    "GA/TC": (1.5, -1.8),
    "GA/TG": (1.9, -1.8),
    "GC/AA": (4.6, 7.5),
    "GC/AC": (4.7, 7.5),
    "GC/AT": (5.2, 7.5),
    "GC/GA": (4.1, 5.9),
    "GC/GC": (4.2, 5.9),
    "GC/GT": (4.6, 5.9),
    "GC/TA": (2.7, -0.6),
    "GC/TC": (2.3, -0.6),
    "GC/TT": (2.1, -0.6),
    "GG/AA": (4.1, 7.0),
    "GG/AG": (4.1, 7.0),
    "GG/AT": (4.7, 7.0),
    "GG/GA": (3.6, 5.4),
    "GG/GG": (3.6, 5.4),
    "GG/GT": (4.1, 5.4),
    "GG/TA": (2.2, -1.1),
    "GG/TG": (2.2, -1.1),
    "GG/TT": (1.6, -1.1),
    "GT/AC": (4.6, 6.6),
    "GT/AG": (4.5, 6.6),
    "GT/AT": (5.1, 6.6),
    "GT/GC": (4.2, 5.0),
    "GT/GG": (4.1, 5.0),
    "GT/GT": (4.6, 5.0),
    "GT/TC": (2.2, -1.5),
    "GT/TG": (2.6, -1.5),
    "GT/TT": (2.0, -1.5),
    "TA/CA": (4.6, 6.6),
    "TA/CC": (4.8, 6.6),
    "TA/CG": (5.1, 6.6),
    "TA/GA": (4.0, 4.7),
    "TA/GC": (4.1, 4.7),
    "TA/GG": (4.0, 4.7),
    "TA/TA": (2.5, -1.8),
    "TA/TC": (2.1, -1.8),
    "TA/TG": (2.5, -1.8),
    "TC/CA": (4.6, 7.8),
    "TC/CC": (4.8, 7.8),
    "TC/CT": (4.7, 7.8),
    "TC/GA": (3.9, 5.9),
    "TC/GC": (4.0, 5.9),
    "TC/GT": (4.4, 5.9),
    "TC/TA": (2.5, -0.6),
    "TC/TC": (2.1, -0.6),
    "TC/TT": (1.9, -0.6),
    "TG/CA": (4.9, 7.3),
    "TG/CG": (5.4, 7.3),
    "TG/CT": (5.0, 7.3),
    "TG/GA": (4.2, 5.4),
    "TG/GG": (4.2, 5.4),
    "TG/GT": (4.7, 5.4),
    "TG/TA": (2.8, -1.1),
    "TG/TG": (2.8, -1.1),
    "TG/TT": (2.2, -1.1),
    "TT/CC": (4.3, 6.9),
    "TT/CG": (4.6, 6.9),
    "TT/CT": (4.2, 6.9),
    "TT/GC": (3.6, 5.0),
    "TT/GG": (3.5, 5.0),
    "TT/GT": (4.0, 5.0),
    "TT/TC": (1.6, -1.5),
    "TT/TG": (2.0, -1.5),
    "TT/TT": (1.4, -1.5),
}


def run_melting_batch(
    sequences,
    complements=None,
    filename="melting_batch_input.txt",
    Na=0.075,
    Mg=0.01,
    Tris=0.02,
    formamide=20,
    K=0.075,
    dnac=0.00000006,
    form_method="bla96",
    oligo_conc_form=4,
    use_slurm=False,
):
    """
    Run melting-batch command for a batch of sequences and return their melting temperatures.

    Parameters:
        sequences (list): List of sequences (strings) to process.
        complements (list): List of complementary sequences (optional). If not provided, they will be generated automatically.
        filename (str): Name of the batch input file to create.
        Na, Mg, Tris, formamide, K, dnac: Experimental concentration parameters for the calculation, applied to all sequences.
        form_method (str): Formamide correction method to use. Choose either 'lincorr' (% formamide) or 'bla96' (mol/L)
        oligo_conc_form (int): Oligo concentration formula (1 for equimolar oligos, 4 for sequence in excess of complement)
        use_slurm (bool): If True, run melting-batch on the current SLURM node (expects `ml melting` already loaded).

    Returns:
        dict: A dictionary where keys are sequence pairs, and values are their respective melting temperatures.
    """
    # Create the melting-batch input file
    if complements is None:
        complements = [str(Seq(s).complement()) for s in sequences]
    else:
        complements = [str(Seq(c).complement()) for c in complements]

    formamide = (formamide * 10 * 1.13) / 45.04

    # Ensure absolute path for input file
    with open(filename, "w") as f:
        for s_arg, c_arg in zip(sequences, complements):
            # Write each sequence and its complement on the same line, separated by a space
            f.write(f"{s_arg} {c_arg}\n")
    # Get the full path to the file
    filename = os.path.abspath(filename)
    filename_base = os.path.basename(filename)

    args = [
        "-for",
        f"{form_method}",
        "-F",
        f"{oligo_conc_form}",
        "-E",
        f"Na={Na}:Mg={Mg}:Tris={Tris}:formamide={formamide}:K={K}",
        "-P",
        str(dnac),
        "-H",
        "dnadna",
        filename,
    ]

    cmd = "melting-batch"

    # Run the melting-batch command (SLURM vs Windows)
    if use_slurm:
        # Run on current node; assume module is already loaded and cmd is on PATH
        completed = subprocess.run(
            [cmd] + args,
            capture_output=True,
            text=True,
            check=False,
        )
        used_cwd = os.getcwd()
    else:
        # Original Windows behavior
        working_dir = "C:/ProgramData/MELTING5.2.0/executable/"
        completed = subprocess.run(
            [cmd] + args,
            cwd=working_dir,
            capture_output=True,
            text=True,
            shell=True,
            check=False,
        )
        used_cwd = working_dir

    if completed.returncode != 0:
        raise RuntimeError(
            f"melting-batch failed (code {completed.returncode}). stderr:\n{completed.stderr}"
        )

    # Locate the result file robustly
    candidates = [
        f"{filename}.results.csv",
        os.path.join(used_cwd, f"{filename_base}.results.csv"),
        os.path.join(os.getcwd(), f"{filename_base}.results.csv"),
    ]
    results_path = next((p for p in candidates if os.path.isfile(p)), None)
    if results_path is None:
        raise FileNotFoundError(
            "Unable to locate melting-batch results CSV. Tried:\n"
            + "\n".join(candidates)
        )

    # Parse the output to extract melting temperatures
    output = pd.read_csv(results_path, sep="\t")

    melting_temps = {}

    # Iterate over the DataFrame rows and extract sequences and melting temperatures
    for _, row in output.iterrows():
        seq_pair = (row["Sequence"], str(Seq(row["Complementary"]).complement()))
        melting_temp = row["Tm (deg C)"]
        if isinstance(melting_temp, str):
            melting_temp = float(melting_temp.replace(",", ""))
        melting_temps[seq_pair] = melting_temp

    return melting_temps


# Modify the process_row function
def process_melting_row(row, armlength=20):
    query_seq = row.qseq
    subject_seq = row.sseq
    ligation_site = armlength + 1
    query_start = row.qstart
    (
        query_left,
        query_right,
        subject_left,
        subject_right,
    ) = split_arms(query_seq, subject_seq, ligation_site, query_start)

    # Check for gaps or mismatches near the ligation site
    if not has_gap_or_mismatch(query_seq, subject_seq, ligation_site, query_start):
        # Split the sequences into left and right arms using query-derived split
        (
            query_left,
            query_right,
            subject_left,
            subject_right,
        ) = split_arms(query_seq, subject_seq, ligation_site, query_start)
        # Fill gaps
        if query_left:
            query_left, subject_left = fill_gaps(query_left, subject_left)
        else:
            query_left, subject_left = None, None
        if query_right:
            query_right, subject_right = fill_gaps(query_right, subject_right)
        else:
            query_right, subject_right = None, None

        return {
            "query_left": query_left,
            "subject_left": subject_left,
            "query_right": query_right,
            "subject_right": subject_right,
            "ligation_site_missmatch": False,
        }
    else:
        return {
            "query_left": query_left,
            "subject_left": subject_left,
            "query_right": query_right,
            "subject_right": subject_right,
            "ligation_site_missmatch": True,
        }


# Batch processing function
def process_dataframe_in_batches(
    df, precomputed_variants, armlength=20, tm_threshold=37, batch_size=300
):
    sequences_left = []
    complements_left = []
    sequences_right = []
    complements_right = []

    batch_results_left = {}
    batch_results_right = {}

    melting_results_NN = []

    # tqdm for showing progress
    for row in tqdm(df.itertuples(), total=df.shape[0], desc="Calculating Tms"):
        processed_row = process_row(
            row, armlength, precomputed_variants=precomputed_variants
        )

        if processed_row:
            query_left = processed_row["query_left"]
            subject_left = processed_row["subject_left"]
            query_right = processed_row["query_right"]
            subject_right = processed_row["subject_right"]
            ligation_site_missmatch = processed_row["ligation_site_missmatch"]

            # Append to batch lists
            sequences_left.append(query_left)
            complements_left.append(subject_left)
            sequences_right.append(query_right)
            complements_right.append(subject_right)

            # Run calc_tm_NN for the left and right sequences
            try:
                tm_left_NN = calc_tm_NN(query_left, subject_left)
            except Exception as e:
                tm_left_NN = None

            try:
                tm_right_NN = calc_tm_NN(query_right, subject_right)
            except Exception as e:
                tm_right_NN = None

            # Add results to the list
            melting_results_NN.append(
                {
                    "tm_left_NN": tm_left_NN,
                    "tm_right_NN": tm_right_NN,
                    "left_length": len(query_left),
                    "right_length": len(query_right),
                    "ligation_site_missmatch": ligation_site_missmatch,
                }
            )

            # Check if we have reached the batch size limit
            if len(sequences_left) == batch_size:
                # Run melting for the current batch
                batch_results_left.update(
                    run_melting_batch(sequences_left, complements_left)
                )
                batch_results_right.update(
                    run_melting_batch(sequences_right, complements_right)
                )

                # Clear the batches
                sequences_left.clear()
                complements_left.clear()
                sequences_right.clear()
                complements_right.clear()
        else:
            # Ensure that melting_results_NN is updated even when processed_row is None
            melting_results_NN.append(
                {
                    "tm_left_NN": None,
                    "tm_right_NN": None,
                    "left_length": 0,
                    "right_length": 0,
                    "ligation_site_missmatch": None,
                }
            )

    # Process any remaining sequences in the batch
    if sequences_left:
        batch_results_left.update(run_melting_batch(sequences_left, complements_left))
        batch_results_right.update(
            run_melting_batch(sequences_right, complements_right)
        )

    # Now, iterate through the DataFrame again to assign melting temperatures
    melting_results_combined = []
    for i, row in enumerate(
        tqdm(df.itertuples(), total=df.shape[0], desc="Collecting results")
    ):
        processed_row = process_melting_row(row, armlength)
        if processed_row:
            query_left = processed_row["query_left"]
            subject_left = processed_row["subject_left"]
            query_right = processed_row["query_right"]
            subject_right = processed_row["subject_right"]
            # Fetch Tm values from batch results
            tm_left = batch_results_left.get((query_left, subject_left), None)
            tm_right = batch_results_right.get((query_right, subject_right), None)
            # convert from strings to floats
            tm_left = float(tm_left) if tm_left is not None else None
            tm_right = float(tm_right) if tm_right is not None else None

            # Check if both Tm values are valid
            if (tm_left is not None) and (tm_right is not None):
                valid_probe = (
                    (tm_left >= tm_threshold)
                    and (tm_right >= tm_threshold)
                    and (row.ligation_site_missmatch == False)
                )
            else:
                valid_probe = False

            tm_left_NN = melting_results_NN[i]["tm_left_NN"]
            tm_right_NN = melting_results_NN[i]["tm_right_NN"]

            if (tm_left_NN is not None) and (tm_right_NN is not None):
                valid_probe_NN = (
                    (tm_left_NN >= tm_threshold)
                    and (tm_right_NN >= tm_threshold)
                    and (row.ligation_site_missmatch == False)
                )
            else:
                valid_probe_NN = False

            melting_results_combined.append(
                {
                    "tm_left_melting": tm_left,
                    "tm_right_melting": tm_right,
                    "valid_probe_melting": valid_probe,
                    "tm_left_NN": tm_left_NN,
                    "tm_right_NN": tm_right_NN,
                    "valid_probe_NN": valid_probe_NN,
                    "ligation_site_missmatch": melting_results_NN[i][
                        "ligation_site_missmatch"
                    ],
                    "left_length": melting_results_NN[i]["left_length"],
                    "right_length": melting_results_NN[i]["right_length"],
                }
            )
        else:
            melting_results_combined.append(
                {
                    "tm_left_melting": None,
                    "tm_right_melting": None,
                    "valid_probe_melting": None,
                    "tm_left_NN": None,
                    "tm_right_NN": None,
                    "valid_probe_NN": None,
                    "ligation_site_missmatch": None,
                    "left_length": 0,
                    "right_length": 0,
                }
            )

    # Convert results to a DataFrame
    results_df = pd.DataFrame(melting_results_combined)

    df = df.reset_index(drop=True)
    results_df = results_df.reset_index(drop=True)

    # Concatenate the results with the original DataFrame
    for col in results_df.columns:
        df.loc[:, col] = results_df[col]

    for i in range(20, 51):
        df[f"valid_{i}"] = (
            (df["tm_left_melting"] > i)
            & (df["tm_right_melting"] > i)
            & (~df["ligation_site_missmatch"])
        )

    return df
