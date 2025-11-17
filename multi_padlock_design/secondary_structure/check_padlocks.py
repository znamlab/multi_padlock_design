from __future__ import annotations

import math
import os
import pickle
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
from collections import defaultdict
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.Seq import Seq

try:
    import multi_padlock_design.config as config
except ModuleNotFoundError:
    from multi_padlock_design import config

from tqdm import tqdm
from tqdm.auto import tqdm as _tqdm
from znamutils import slurm_it

from multi_padlock_design.blast.readblast import (
    calc_tm_NN,
    fill_gaps,
    has_gap_or_mismatch,
    split_arms,
)
from multi_padlock_design.io import formatrefseq


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
    # Trigger variant discovery (result unused here; retained for side effects or future use)
    find_variants(gene)
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
    "GT/AC": (4.6, 6.6),
    "GT/AG": (4.5, 6.6),
    "GT/AT": (5.1, 6.6),
    "GT/GC": (4.2, 5.0),
    "GT/GG": (4.1, 5.0),
    "GT/GT": (4.6, 5.0),
    "GT/TC": (2.2, -1.5),
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
            except Exception:
                tm_left_NN = None

            try:
                tm_right_NN = calc_tm_NN(query_right, subject_right)
            except Exception:
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
                    and (not row.ligation_site_missmatch)
                )
            else:
                valid_probe = False

            tm_left_NN = melting_results_NN[i]["tm_left_NN"]
            tm_right_NN = melting_results_NN[i]["tm_right_NN"]

            if (tm_left_NN is not None) and (tm_right_NN is not None):
                valid_probe_NN = (
                    (tm_left_NN >= tm_threshold)
                    and (tm_right_NN >= tm_threshold)
                    and (not row.ligation_site_missmatch)
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


# -----------------------------
# User-configurable hetero/homo/hairpin parameters
# -----------------------------
KMER = 10
MAX_CANDIDATES = 250
MAX_FOR_THERMO = 80
ALLOW_SELF_AS_PARTNER = False

TEMP_C = 37  # temperature in °C for RNAstructure
FOLD_BIN = os.environ.get("RNASTRUCTURE_FOLD_BIN", "Fold")
DUPLEX_BIN = os.environ.get("RNASTRUCTURE_DUPLEX_BIN", "DuplexFold")

# DATAPATH for RNAstructure data tables
RNASTRUCTURE_DATAPATH = os.environ.get(
    "DATAPATH",
    os.path.join(
        os.environ.get("CONDA_PREFIX", ""), "share", "rnastructure", "data_tables"
    ),
)


def _env_with_datapath():
    """
    Get the environment variables with the RNAstructure DATAPATH.
    """
    env = os.environ.copy()
    if RNASTRUCTURE_DATAPATH:
        env["DATAPATH"] = RNASTRUCTURE_DATAPATH
    return env


# -----------------------------
# Utilities
# -----------------------------
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    """
    Get the reverse complement of a DNA sequence.

    Args:
        seq (str): The DNA sequence to reverse complement.

    Returns:
        str: The reverse complement of the input DNA sequence.
    """
    return seq.translate(_COMP)[::-1]


def clean_seq(seq: str) -> str:
    """
    Clean the input DNA sequence by removing whitespace and converting to uppercase.

    Args:
        seq (str): The DNA sequence to clean.

    Returns:
        str: The cleaned DNA sequence.
    """
    s = (seq or "").strip().upper()
    s = "".join(c if c in "ACGT" else "N" for c in s)
    return s


# Memoized per-process binary check
_BINS_CHECKED = False


def ensure_bins():
    """Check RNAstructure binaries once per process.

    Raises:
        RuntimeError: If required binaries are not found.
    """
    global _BINS_CHECKED
    if _BINS_CHECKED:
        return
    for b in (FOLD_BIN, DUPLEX_BIN):
        try:
            subprocess.run(
                [b, "--help"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            )
        except FileNotFoundError:
            raise RuntimeError(
                f"Required RNAstructure binary '{b}' not found on PATH. "
                "Install via 'conda install rnastructure' or set RNASTRUCTURE_*_BIN env vars."
            )
    _BINS_CHECKED = True


def write_fasta(seq: str, path: str, name="seq"):
    """Write a FASTA file.

    Args:
        seq (str): The DNA sequence to write.
        path (str): The file path to write the FASTA file.
        name (str, optional): The name for the FASTA entry. Defaults to "seq".
    """
    with open(path, "w") as f:
        f.write(f">{name}\n{seq}\n")


def parse_ct_energy(ct_path: str):
    """
    Parse RNAstructure CT header and return energy (kcal/mol, negative is favorable).

    Args:
        ct_path (str): Path to CT file.
    Returns:
        float: Energy in kcal/mol, or np.nan if not found.
    """
    try:
        with open(ct_path, "r") as f:
            header = f.readline()
        m = re.search(
            r"(?:ENERGY|dG)\s*=\s*([+-]?\d+(?:\.\d+)?)", header, flags=re.IGNORECASE
        )
        if m:
            return float(m.group(1))
    except Exception:
        pass
    return np.nan


def rnastructure_hairpin_dg(seq: str, temp_c: float = TEMP_C) -> float:
    """
    Single-strand MFE (hairpins included) via RNAstructure 'Fold' with --DNA.
    Returns kcal/mol

    Args:
        seq (str): DNA sequence.
        temp_c (float): Temperature in Celsius (default: 37).
    Returns:
        float: ΔG of folding in kcal/mol.
    """
    ensure_bins()
    tmpdir = tempfile.mkdtemp(prefix="rnastruct_hp_")
    try:
        fa = os.path.join(tmpdir, "s.fa")
        ct = os.path.join(tmpdir, "out.ct")
        write_fasta(seq, fa)
        # Temperature in Kelvin for RNAstructure
        temp_k = 273.15 + temp_c
        cmd = [FOLD_BIN, fa, ct, "--DNA", "--MFE", "--temperature", f"{temp_k}"]
        res = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=_env_with_datapath(),
        )
        if res.returncode != 0:
            raise RuntimeError(
                f"Fold failed: {res.stderr.strip() or res.stdout.strip()}"
            )
        return parse_ct_energy(ct)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def rnastructure_duplex_dg(seq1: str, seq2: str, temp_c: float = TEMP_C) -> float:
    """
    Intermolecular duplex MFE via RNAstructure 'DuplexFold' with --DNA.
    Returns kcal/mol. No intramolecular pairs are allowed (ideal for dimers).

    Args:
        seq1 (str): First DNA sequence.
        seq2 (str): Second DNA sequence.
        temp_c (float): Temperature in Celsius (default: 37).
    Returns:
        float: ΔG of duplex formation in kcal/mol.
    """
    ensure_bins()
    tmpdir = tempfile.mkdtemp(prefix="rnastruct_duplex_")
    try:
        fa1 = os.path.join(tmpdir, "a.fa")
        fa2 = os.path.join(tmpdir, "b.fa")
        ct = os.path.join(tmpdir, "out.ct")
        write_fasta(seq1, fa1, "a")
        write_fasta(seq2, fa2, "b")
        temp_k = 273.15 + temp_c
        cmd = [DUPLEX_BIN, fa1, fa2, ct, "--DNA", "--temperature", f"{temp_k}"]
        res = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=_env_with_datapath(),
        )
        if res.returncode != 0:
            raise RuntimeError(
                f"DuplexFold failed: {res.stderr.strip() or res.stdout.strip()}"
            )
        return parse_ct_energy(ct)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def longest_exact_rc_run(a: str, b: str) -> int:
    """
    Finds longest exact reverse-complement match between a and b.

    Args:
        a (str): First sequence.
        b (str): Second sequence.
    Returns:
        int: Length of longest exact reverse-complement match.
    """
    rb = revcomp(b)
    n, m = len(a), len(rb)
    prev = [0] * (m + 1)
    best = 0
    for i in range(1, n + 1):
        curr = [0] * (m + 1)
        ai = a[i - 1]
        for j in range(1, m + 1):
            if ai == rb[j - 1]:
                curr[j] = prev[j - 1] + 1
                if curr[j] > best:
                    best = curr[j]
        prev = curr
    return best


def three_prime_bias_weight(pos: int, L: int, k: int) -> float:
    """Calculate the 3' bias weight for a given position in a sequence.

    Args:
        pos (int): The position of interest (0-based).
        L (int): The total length of the sequence.
        k (int): The length of the k-mer.

    Returns:
        float: The 3' bias weight.
    """
    dist3 = L - (pos + k)
    return 1.0 / (1.0 + dist3)


def _intra_task(item):
    """Compute hairpin and homodimer ΔG values for a given sequence.

    Args:
        item (tuple): A tuple containing the index and the sequence.

    Returns:
        tuple: A tuple containing the index, hairpin ΔG, and homodimer Δ
    """
    i, seq = item
    hp = np.nan
    hd = np.nan
    try:
        hp = rnastructure_hairpin_dg(seq)
    except Exception:
        pass
    try:
        hd = rnastructure_duplex_dg(seq, seq)
    except Exception:
        pass
    return i, hp, hd


# -----------------------------
# Build / normalize input df
# -----------------------------
def prepare_df(df: pd.DataFrame) -> pd.DataFrame:
    """Prepare the DataFrame by adding a name column and cleaning sequences.

    Args:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        pd.DataFrame: The prepared DataFrame with cleaned sequences and a name column.
    """
    if "name" not in df.columns:
        df = df.copy()
        df["name"] = [f"oligo_{i}" for i in range(len(df))]
    df = df.copy()
    df["padlock"] = df["padlock"].map(clean_seq)
    bad = df["padlock"].str.contains("N")
    if bad.any():
        print(f"[warn] Dropping {bad.sum()} sequences containing non-ACGT symbols.")
        df = df.loc[~bad].reset_index(drop=True)
    return df


# -----------------------------
# Hairpin & homodimer (RNAstructure)
# -----------------------------
def compute_intrastrand(df: pd.DataFrame, n_jobs: int | None = None) -> pd.DataFrame:
    """
    Compute hairpin and homodimer ΔG values with accurate, single-line tqdm updates.

    Args:
        df (pd.DataFrame): Input DataFrame with a "padlock" column containing sequences
        n_jobs (int | None): Number of parallel jobs to use. If None or 1, runs serially.
    Returns:
        pd.DataFrame: DataFrame with added "hairpin_dg_kcalmol" and
                      "homodimer_dg_kcalmol" columns.
    """
    seqs = df["padlock"].tolist()
    N = len(seqs)
    hairpin = [np.nan] * N
    homodimer = [np.nan] * N

    if n_jobs is None or n_jobs == 1:
        pbar = _pbar(total=N, desc="Hairpin & homodimer", unit="oligo")
        for i, s in enumerate(seqs):
            _, hp, hd = _intra_task((i, s))
            hairpin[i] = hp
            homodimer[i] = hd
            pbar.update(1)
        pbar.close()
    else:
        with ProcessPoolExecutor(max_workers=n_jobs, initializer=ensure_bins) as ex:
            pbar = _pbar(total=N, desc="Hairpin & homodimer (parallel)", unit="oligo")
            it = iter(enumerate(seqs))
            inflight = set()
            max_inflight = n_jobs * 4

            for _ in range(min(N, max_inflight)):
                try:
                    inflight.add(ex.submit(_intra_task, next(it)))
                except StopIteration:
                    break

            while inflight:
                done, inflight = wait(inflight, return_when=FIRST_COMPLETED)
                for fut in done:
                    i, hp, hd = fut.result()
                    hairpin[i] = hp
                    homodimer[i] = hd
                    pbar.update(1)
                while len(inflight) < max_inflight:
                    try:
                        inflight.add(ex.submit(_intra_task, next(it)))
                    except StopIteration:
                        break
            pbar.close()

    out = df.copy()
    out["hairpin_dg_kcalmol"] = hairpin
    out["homodimer_dg_kcalmol"] = homodimer
    return out


# -----------------------------
# k-mer index on forward seqs
# -----------------------------
def build_kmer_index(seqs, k: int):
    """Build a k-mer index from the given sequences.

    Args:
        seqs (list): A list of sequences.
        k (int): The length of the k-mer.

    Returns:
        dict: A dictionary mapping k-mers to their positions in the sequences.
    """
    index = defaultdict(list)  # kmer -> [(idx, pos)]
    for idx, s in enumerate(seqs):
        L = len(s)
        for pos in range(0, L - k + 1):
            index[s[pos : pos + k]].append((idx, pos))
    return index


def generate_candidates_for_oligo(i, seqs, index, k: int, max_candidates: int):
    """Generate candidate oligos for a given oligo index.

    Args:
        i (int): The index of the oligo to generate candidates for.
        seqs (list): The list of oligo sequences.
        index (dict): The k-mer index.
        k (int): The k-mer length.
        max_candidates (int): The maximum number of candidate oligos to return.

    Returns:
        list: A list of candidate oligo indices.
    """
    s = seqs[i]
    L = len(s)
    candidate_scores = defaultdict(float)
    for pos_i in range(0, L - k + 1):
        rc = revcomp(s[pos_i : pos_i + k])
        hits = index.get(rc)
        if not hits:
            continue
        w_i = three_prime_bias_weight(pos_i, L, k)
        for j, pos_j in hits:
            if (not ALLOW_SELF_AS_PARTNER) and j == i:
                continue
            Lj = len(seqs[j])
            w_j = three_prime_bias_weight(pos_j, Lj, k)
            combined = 0.6 * w_i + 0.4 * w_j
            if combined > candidate_scores[j]:
                candidate_scores[j] = combined
    if not candidate_scores:
        return []
    top = sorted(candidate_scores.items(), key=lambda kv: kv[1], reverse=True)[
        :max_candidates
    ]
    return [j for j, _ in top]


def prescore_pairs_by_runlen(i, candidate_js, seqs, max_for_thermo: int):
    """Prescore candidate pairs by longest exact reverse-complement run length.
    Args:
        i (int): The index of the oligo to prescore candidates for.
        candidate_js (list): The list of candidate oligo indices.
        seqs (list): The list of oligo sequences.
        max_for_thermo (int): The maximum number of candidates to return for thermodynamic evaluation.
    Returns:
        list: A list of prescored candidate oligo indices.
    """
    if not candidate_js:
        return []
    s_i = seqs[i]
    scored = []
    for j in candidate_js:
        Lrun = longest_exact_rc_run(s_i, seqs[j])
        scored.append((j, Lrun, 1000 * Lrun + (j % 7) * 0.01))
    scored.sort(key=lambda x: x[2], reverse=True)
    return [j for j, _, _ in scored[:max_for_thermo]]


# -----------------------------
# Parallel pool: worker state
# -----------------------------
_G_SEQS = None
_G_NAMES = None


def _init_pool(seqs, names):
    """Initializer run in each worker process to stash read-only data & warm binaries.

    Args:
        seqs (list): List of sequences.
        names (list): List of names.
    """
    global _G_SEQS, _G_NAMES
    _G_SEQS = seqs
    _G_NAMES = names
    ensure_bins()


def _duplex_task(pair: Tuple[int, int]):
    """Compute ΔG for a single unordered pair (i, j).

    Args:
        pair (tuple): A tuple containing two indices (i, j).
    Returns:
        tuple: A tuple containing the pair (i, j) and the computed ΔG value
    """
    i, j = pair
    try:
        dg = rnastructure_duplex_dg(_G_SEQS[i], _G_SEQS[j])
    except Exception:
        dg = np.nan
    return (i, j), dg


def _pbar(total: int, desc: str, unit: str):
    """Create a stable, single-line progress bar (or a no-op in non-TTY).

    Args:
        total (int): Total number of iterations.
        desc (str): Description for the progress bar.
        unit (str): Unit of measurement for the progress bar.
    Returns:
        tqdm: A tqdm progress bar instance.
    """
    return _tqdm(
        total=total,
        desc=desc,
        unit=unit,
        dynamic_ncols=True,
        mininterval=0.2,
        smoothing=0.1,
        leave=False,  # <- do not leave bar lines behind
        position=0,  # <- keep on a single line
    )


# -----------------------------
# End-to-end driver (parallel + dedup + tqdm)
# -----------------------------
def annotate_with_thermo(
    df: pd.DataFrame,
    kmer=KMER,
    max_candidates=MAX_CANDIDATES,
    max_for_thermo=MAX_FOR_THERMO,
    n_jobs: int | None = None,  # control parallelism
) -> pd.DataFrame:
    """
    Annotate with hairpin/homodimer (serial) and best heterodimer (parallelized with dedup).
    Each unique unordered pair (i, j) from all shortlists is evaluated exactly once.

    Args:
        df (pd.DataFrame): Input DataFrame with a "padlock" column containing sequences
        kmer (int): k-mer length for indexing (default: 10)
        max_candidates (int): Max candidates to shortlist per oligo (default: 250)
        max_for_thermo (int): Max candidates per oligo to evaluate thermodynamically (default: 80)
        n_jobs (int | None): Number of parallel jobs to use. If None, defaults to min(8, cpu_count).
    Returns:
        pd.DataFrame: DataFrame with added columns:
                      "hairpin_dg_kcalmol", "homodimer_dg_kcalmol",
                      "best_heterodimer_dg_kcalmol", "best_heterodimer_partner",
                      "heterodimer_candidate_dgs_kcalmol", "heterodimer_candidate_partners"
    """
    df = prepare_df(df)

    # Step 1: Hairpin & homodimer (serial)
    df1 = compute_intrastrand(df, n_jobs=n_jobs)

    names = df1["name"].tolist()
    seqs = df1["padlock"].tolist()
    N = len(seqs)

    # Step 2: k-mer index
    kindex = build_kmer_index(seqs, kmer)

    # Step 3a: per-oligo candidate shortlist (serial, fast; show progress)
    shortlists: List[List[int]] = [None] * N
    pbar = _pbar(total=N, desc="Shortlisting partners", unit="oligo")
    for i in range(N):
        cand_js = generate_candidates_for_oligo(i, seqs, kindex, kmer, max_candidates)
        shortlist_js = prescore_pairs_by_runlen(i, cand_js, seqs, max_for_thermo)
        shortlists[i] = shortlist_js
        pbar.update(1)
    pbar.close()

    # Step 3b: build set of unique unordered pairs across ALL shortlists
    # (ΔG(i,j) == ΔG(j,i)); evaluate each pair exactly once.
    seen = set()
    pairs: List[Tuple[int, int]] = []
    for i, js in enumerate(shortlists):
        if not js:
            continue
        for j in js:
            a, b = (i, j) if i < j else (j, i)
            if (not ALLOW_SELF_AS_PARTNER) and a == b:
                continue
            if (a, b) not in seen:
                seen.add((a, b))
                pairs.append((a, b))

    # Early exit: no pairs to evaluate
    if not pairs:
        out = df1.copy()
        out["best_heterodimer_dg_kcalmol"] = [np.nan] * N
        out["best_heterodimer_partner"] = [None] * N
        return out

    # Step 3c: parallel DuplexFold over unique pairs with tqdm and bounded in-flight futures
    if n_jobs is None:
        n_jobs = max(1, min(8, (os.cpu_count() or 2)))

    energy: Dict[Tuple[int, int], float] = {}
    with ProcessPoolExecutor(
        max_workers=n_jobs, initializer=_init_pool, initargs=(seqs, names)
    ) as ex:
        pbar = _pbar(total=len(pairs), desc="DuplexFold (unique pairs)", unit="pair")
        inflight = set()
        it = iter(pairs)
        max_inflight = n_jobs * 4

        # Prime submissions
        for _ in range(min(len(pairs), max_inflight)):
            try:
                pair = next(it)
            except StopIteration:
                break
            inflight.add(ex.submit(_duplex_task, pair))

        # Consume as they finish; keep queue topped up
        while inflight:
            done, inflight = wait(inflight, return_when=FIRST_COMPLETED)
            for fut in done:
                (i, j), dg = fut.result()
                energy[(i, j)] = dg
                pbar.update(1)
            while len(inflight) < max_inflight:
                try:
                    pair = next(it)
                except StopIteration:
                    break
                inflight.add(ex.submit(_duplex_task, pair))
        pbar.close()

    # Step 3d: for each oligo, pick best partner from its shortlist using cached energies
    best_hdgs = [np.nan] * N
    best_partners = [None] * N
    cand_dgs = [[] for _ in range(N)]
    cand_names = [[] for _ in range(N)]

    pbar = _pbar(total=N, desc="Ranking partners", unit="oligo")
    for i, js in enumerate(shortlists):
        ranked: list[tuple[float, int]] = []
        if js:
            for j in js:
                a, b = (i, j) if i < j else (j, i)
                dg = energy.get((a, b), np.nan)
                if not np.isnan(dg):
                    ranked.append((float(dg), j))

            # Sort strongest (most negative ΔG) -> weakest
            ranked.sort(key=lambda x: x[0])

            # Save full ranked lists
            cand_dgs[i] = [dg for dg, _j in ranked]
            cand_names[i] = [names[_j] for _dg, _j in ranked]

            # Also keep the legacy "best_*" summaries for convenience/back-compat
            if ranked:
                best_hdgs[i] = ranked[0][0]
                best_partners[i] = names[ranked[0][1]]

        pbar.update(1)
    pbar.close()

    out = df1.copy()
    out["best_heterodimer_dg_kcalmol"] = best_hdgs
    out["best_heterodimer_partner"] = best_partners
    # NEW: full ranked candidate lists per oligo
    out["heterodimer_candidate_dgs_kcalmol"] = cand_dgs
    out["heterodimer_candidate_partners"] = cand_names

    return out


# =============================================================
# Exhaustive heterodimer thermodynamics (all unique unordered pairs)
# Batched across SLURM jobs using @slurm_it
# =============================================================


@slurm_it(
    conda_env="multi-padlock-design",
    from_imports={"lib.check_padlocks": ("run_exhaustive_heterodimer_batch")},
)
def run_exhaustive_heterodimer_batch(
    probe_df_path: str,
    output_dir: str,
    batch_index: int,
    n_batches: int = 1000,
    sequence_col: str = "padlock",
    name_col: Optional[str] = None,
    overwrite: bool = False,
):
    """Compute duplex ΔG for one batch of an exhaustive enumeration of all
    unique unordered probe pairs.

    Pair indexing strategy:
        For N sequences, total M = N*(N-1)/2 pairs (i<j) are linearized
        in lexicographic order of (i,j). The interval [0, M) is split
        into n_batches ~equal contiguous chunks; this batch processes
        the slice assigned by batch_index.

    Output:
        A pickle containing a DataFrame with columns:
            i, j, probe_i, probe_j, dg
        saved to: {output_dir}/heterodimer_pairs_part_{batch_index:04d}.pkl
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_file = output_dir / f"heterodimer_pairs_part_{batch_index:04d}.pkl"
    if out_file.exists() and not overwrite:
        return str(out_file)

    # Load probe dataframe (csv or pickle)
    probe_df_path = Path(probe_df_path)
    if probe_df_path.suffix.lower() in {".pkl", ".pickle"}:
        df = pd.read_pickle(probe_df_path)
    else:
        df = pd.read_csv(probe_df_path)

    if sequence_col not in df.columns:
        raise ValueError(f"Sequence column '{sequence_col}' not found")

    if name_col is None:
        # try common name columns
        for cand in ["name", "probe_name", "id", "gene"]:
            if cand in df.columns:
                name_col = cand
                break
    if name_col is None:
        name_col = "_autoname"
        df = df.copy()
        df[name_col] = [f"probe_{i}" for i in range(len(df))]

    seqs = df[sequence_col].map(clean_seq).tolist()
    names = df[name_col].tolist()
    N = len(seqs)
    if N < 2:
        raise ValueError("Need at least two probes")

    total_pairs = N * (N - 1) // 2
    n_batches = min(n_batches, total_pairs) if total_pairs else 0
    if n_batches == 0:
        raise ValueError("No pairs to enumerate")
    if batch_index >= n_batches:
        # Over-provisioned batch index: write empty file sentinel
        with out_file.open("wb") as f:
            pickle.dump(pd.DataFrame([]), f)
        return str(out_file)

    # Determine k-range for this batch
    per_batch = math.ceil(total_pairs / n_batches)
    start_k = batch_index * per_batch
    end_k = min(start_k + per_batch, total_pairs)
    if start_k >= end_k:
        with out_file.open("wb") as f:
            pickle.dump(pd.DataFrame([]), f)
        return str(out_file)

    try:
        ensure_bins()
    except Exception:
        pass  # allow fallback heuristic if binaries unavailable

    def pair_generator(N: int, start_k: int, end_k: int):
        """Yield (i,j) pairs for linear indices in [start_k,end_k)."""
        # Advance to starting k
        k = 0
        i = 0
        remaining = start_k
        while i < N - 1:
            block = N - i - 1
            if remaining < block:
                j = i + 1 + remaining
                k = start_k
                break
            remaining -= block
            i += 1
        else:
            return
        while k < end_k and i < N - 1:
            yield i, j
            k += 1
            j += 1
            if j >= N:
                i += 1
                j = i + 1

    # ---------------------------------------------------------
    # Progress over this batch's share of pairs
    # ---------------------------------------------------------
    batch_total = end_k - start_k
    # For SLURM .out logs we avoid overly chatty output: set a min interval.
    # tqdm will rewrite the same line on TTY; on non‑TTY (captured) it still
    # prints periodic refreshes controlled by mininterval.
    # If the stream is not a TTY we still keep it enabled so the .out file
    # accumulates occasional progress lines.
    try:
        _ = sys.stdout.isatty()
    except Exception:
        pass
    # Heuristic miniters: aim for at most ~200 updates.
    miniters = max(1, batch_total // 200) if batch_total > 0 else 1
    pbar = tqdm(
        total=batch_total,
        desc=f"heterodimer batch {batch_index+1}/{n_batches}",
        unit="pair",
        smoothing=0.0,
        mininterval=30.0,  # seconds between forced refreshes
        miniters=miniters,
        dynamic_ncols=False,
        disable=False,  # keep enabled to write progress lines to SLURM log
        ascii=True,
        leave=False,
    )

    records = []
    for i, j in pair_generator(N, start_k, end_k):
        s1 = seqs[i]
        s2 = seqs[j]

        dg = rnastructure_duplex_dg(s1, s2)

        records.append(
            {"i": i, "j": j, "probe_i": names[i], "probe_j": names[j], "dg": dg}
        )
        pbar.update(1)
    pbar.close()

    part_df = pd.DataFrame(records)
    with out_file.open("wb") as f:
        pickle.dump(part_df, f)
    return str(out_file)


@slurm_it(
    conda_env="multi-padlock-design",
    from_imports={"lib.check_padlocks": ("aggregate_exhaustive_heterodimer_results",)},
)
def aggregate_exhaustive_heterodimer_results(
    output_dir: str,
    probe_df_path: str,
    sequence_col: str = "padlock",
    matrix_pickle: Optional[str] = None,
    overwrite: bool = True,
):
    """Aggregate batch part files into a symmetric ΔG matrix and save pickle.

    The pickle contains a dict with keys: matrix (NxN float64),
    names (list), sequence_col.
    """
    output_dir = Path(output_dir)
    if matrix_pickle is None:
        matrix_pickle = output_dir / "heterodimer_dg_matrix.pkl"
    else:
        matrix_pickle = Path(matrix_pickle)
    if matrix_pickle.exists() and not overwrite:
        return str(matrix_pickle)

    probe_df_path = Path(probe_df_path)
    if probe_df_path.suffix.lower() in {".pkl", ".pickle"}:
        df = pd.read_pickle(probe_df_path)
    else:
        df = pd.read_csv(probe_df_path)
    if sequence_col not in df.columns:
        raise ValueError(f"Sequence column '{sequence_col}' not found")
    # Determine name column
    name_col = None
    for c in ["name", "probe_name", "id", "gene", "_autoname"]:
        if c in df.columns:
            name_col = c
            break
    if name_col is None:
        name_col = "_autoname"
        df[name_col] = [f"probe_{i}" for i in range(len(df))]

    names = df[name_col].tolist()
    N = len(names)
    matrix = np.full((N, N), np.nan, dtype=float)

    part_files = sorted(output_dir.glob("heterodimer_pairs_part_*.pkl"))
    if not part_files:
        raise FileNotFoundError("No part files found to aggregate")

    # First pass: compute total rows (for progress total)
    total_rows = 0
    for pf in part_files:
        try:
            with pf.open("rb") as f:
                part_df = pickle.load(f)
        except Exception:
            continue
        if isinstance(part_df, pd.DataFrame):
            total_rows += len(part_df)

    # Progress bar for aggregation; sparse refresh for SLURM logs
    try:
        _ = sys.stdout.isatty()
    except Exception:
        pass
    pbar = tqdm(
        total=total_rows,
        desc="aggregate heterodimer",
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
        it = part_df.itertuples()
        for row in it:
            i = row.i
            j = row.j
            dg = row.dg
            if 0 <= i < N and 0 <= j < N:
                matrix[i, j] = dg
                matrix[j, i] = dg
            pbar.update(1)
    pbar.close()
    np.fill_diagonal(matrix, 0.0)

    payload = {"matrix": matrix, "names": names, "sequence_col": sequence_col}
    with open(matrix_pickle, "wb") as f:
        pickle.dump(payload, f)
    return str(matrix_pickle)


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
    dry_run: bool = False,
):
    """Submit 1000 (or fewer if total pairs < 1000) exhaustive heterodimer jobs
    plus a final aggregation job dependent on all of them.

    Returns (batch_job_ids, aggregate_job_id). If dry_run=True, returns ([], None).
    """
    probe_df_path = Path(probe_df_path)
    if probe_df_path.suffix.lower() in {".pkl", ".pickle"}:
        df = pd.read_pickle(probe_df_path)
    else:
        df = pd.read_csv(probe_df_path)
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

    slurm_opts_batch = {
        "time": time,
        "mem": mem,
        "cpus-per-task": cpus_per_task,
    }
    if partition:
        slurm_opts_batch["partition"] = partition

    slurm_opts_agg = {
        "time": dependency_aggregate_time,
        "mem": dependency_aggregate_mem,
        "cpus-per-task": 1,
    }
    if partition:
        slurm_opts_agg["partition"] = partition

    # Submit batched jobs (one per batch_index) explicitly
    batch_job_ids = []
    for i in range(n_batches):
        jid = run_exhaustive_heterodimer_batch(
            probe_df_path=str(probe_df_path),
            output_dir=str(output_dir),
            batch_index=i,
            n_batches=n_batches,
            sequence_col=sequence_col,
            slurm_folder=slurm_folder,
            scripts_name=f"exhaustive_heterodimer_batch_{i:04d}",
            slurm_options=slurm_opts_batch,
            use_slurm=True,
        )
        batch_job_ids.append(jid)

    # Aggregation job depends on all batch jobs
    aggregate_job_id = aggregate_exhaustive_heterodimer_results(
        output_dir=str(output_dir),
        probe_df_path=str(probe_df_path),
        sequence_col=sequence_col,
        slurm_folder=slurm_folder,
        scripts_name="exhaustive_heterodimer_aggregate",
        job_dependency=batch_job_ids,
        slurm_options=slurm_opts_agg,
        use_slurm=True,
    )
    return batch_job_ids, aggregate_job_id


def main_exhaustive():  # CLI helper
    import argparse

    parser = argparse.ArgumentParser(
        description="Submit exhaustive heterodimer thermodynamics jobs"
    )
    parser.add_argument(
        "probe_df_path", type=str, help="Path to probe dataframe (csv/pkl)"
    )
    parser.add_argument(
        "--output-dir", required=True, type=str, help="Directory for batch outputs"
    )
    parser.add_argument(
        "--slurm-folder",
        required=True,
        type=str,
        help="Existing folder for SLURM scripts & logs",
    )
    parser.add_argument(
        "--n-batches",
        type=int,
        default=1000,
        help="Target number of batches (default 1000)",
    )
    parser.add_argument("--sequence-col", type=str, default="padlock")
    parser.add_argument("--partition", type=str, default=None)
    parser.add_argument("--time", type=str, default="12:00:00")
    parser.add_argument("--mem", type=str, default="8G")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    submit_exhaustive_heterodimer_jobs(
        probe_df_path=args.probe_df_path,
        output_dir=args.output_dir,
        slurm_folder=args.slurm_folder,
        n_batches=args.n_batches,
        time=args.time,
        mem=args.mem,
        partition=args.partition,
        sequence_col=args.sequence_col,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main_exhaustive()
