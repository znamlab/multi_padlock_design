# for multi_padlock_design
# Xiaoyan, 2018

import random
from lib.screenseq import chopseq
from Bio import SeqIO
from difflib import SequenceMatcher
import config
import pandas as pd
import ast
import os
import re
import glob


def correctpos(basepos, targets, targetpos, notMapped, mapTmlist, Tm, siteChopped):
    """Correct fragment coordinates to full length mRNA coordinates"""
    targetposnew = []
    Tmnew = []
    notMappednew = []
    targetsnew = []
    c = -1
    for base in basepos:
        if isinstance(base[0], int):  # only one variant
            c += 1
            targetsnew.append(targets[c])
            targetposnew.append(targetpos[c])
            notMappednew.append(notMapped[c])
            subTm = [tm for i, tm in enumerate(Tm[c]) if i in siteChopped[c]]
            Tmnew.append([subTm[i] for i in mapTmlist[c]])
        else:
            temptargets = []
            temppos = []
            tempnomap = []
            temptm = []

            for subbase in base:
                c += 1
                temptargets = temptargets + targets[c]
                for i, pos in enumerate(targetpos[c]):
                    temppos.append(targetpos[c][i] + subbase[0])
                for i, pos in enumerate(notMapped[c]):
                    tempnomap.append(notMapped[c][i] + subbase[0])
                subTm = [tm for i, tm in enumerate(Tm[c]) if i in siteChopped[c]]
                temptm = temptm + [subTm[i] for i in mapTmlist[c]]

            targetsnew.append(temptargets)
            targetposnew.append(temppos)
            notMappednew.append(tempnomap)
            Tmnew.append(temptm)

        print(f'Not mapped after position correction: {notMappednew}')

    return targetsnew, targetposnew, notMappednew, Tmnew


def assembleprobes(targets, genepars, armlength):
    """Fill backbone sequences"""
    linkers = genepars[1]
    Padlocks = []
    for c, probes in enumerate(targets):
        try:
            linker1 = linkers[c][0]
            linker1[0]
        except:
            linker1 = "ATCGTCGGACTGTAGAACTCTGAACCTGTCG"

        try:
            barcode = linkers[c][1]
            barcode[0]
        except:
            barcode = "NNNNNNNNNN"

        try:
            linker2 = linkers[c][2]
            linker2[0]
        except:
            linker2 = "CT"

        padlocks = []
        for probe in probes:
            padlocks.append(
                probe[armlength:] + linker1 + barcode + linker2 + probe[0:armlength]
            )

        Padlocks.append(padlocks)
    return Padlocks


def removeunmapped(notmapped, targetpos, headers, targets, Tm, probes):
    for i, header in enumerate(headers):
        if len(notmapped[i]):
            for j in list(reversed(range(len(targetpos[i])))):
                if targetpos[i][j] in notmapped[i]:
                    del targets[i][j]
                    del Tm[i][j]
                    del targetpos[i][j]
                    del probes[i][j]
    return (probes, Tm, targetpos, targets)


def selectprobes(n, finals, headers, armlength, outpars):
    """Prioritize probes with no homopolymer sequences and choose randomly n candidates per region

    Args:
        n (int): number of probes per region
        finals (list): list of final probes
        headers (list): list of headers
        armlength (int): arm length
        outpars (list): list of output parameters

    Returns:
        probes (list): list of final probes
        Tm (list): list of Tm values
        targetpos (list): list of target positions
        targets (list): list of targets
        filtered_regions (list): list of binding regions
    """
    input_type = outpars[2]
    tempdir = outpars[1]

    probes = finals[0]
    Tm = finals[1]
    targetpos = finals[2]
    targets = finals[3]
    filtered_regions = [0] * len(headers)
    print(f"Probes: {probes}")
    print(f"Headers: {headers}")
    print(f"Tm Values: {Tm}")
    print(f"Target Positions: {targetpos}")
    print(f"Targets: {targets}")
    print(f"Filtered Regions: {filtered_regions}")

    for i, header in enumerate(headers):
        fullheader = header #keep full header for gene ID
        first = header.split()[0]  # take the first token
        header = first.lstrip(">")  # only strip '>' if it exists
        # Compute probe binding end positions for classification
        target_end = [target + (1 + (2 * armlength)) for target in targetpos[i]]
        regions = []

        # use config.annotation_file (cds/utr3) instead of Ensembl cds/cdna files
        if config.annotation_file:
            anno_df = pd.read_csv(config.annotation_file)
            row = None
            if anno_df is not None and len(targetpos[i]) > 0:
                print("Annotation DataFrame found.")
                # Prefer exact match on gene_symbol column; fall back to gene_name/name if needed
                subset = None
                print(f"Header: {header}")
                if "gene_symbol" in anno_df.columns:
                    subset = anno_df[anno_df["gene_symbol"] == header]
                if (subset is None or subset.empty) and "gene_name" in anno_df.columns:
                    subset = anno_df[anno_df["gene_name"] == header]
                if (subset is None or subset.empty) and "name" in anno_df.columns:
                    subset = anno_df[anno_df["name"] == header]
                # If no symbol available or still no match, leave row as None
                if subset is not None and not subset.empty:
                    row = subset.iloc[0]
            else:
                print("No Annotation DataFrame found.")
            if row is None:
                print("row is None")
                # Cannot classify without a match
                regions = [None] * len(targetpos[i])
            else:
                # Parse CDS tuple "(start, end)"; 3'UTR available but not needed for classification
                cds_tuple = None
                try:
                    val = row["cds"]
                    if pd.notna(val):
                        cds_tuple = ast.literal_eval(str(val))
                except Exception:
                    cds_tuple = None

                if isinstance(cds_tuple, (list, tuple)) and len(cds_tuple) >= 2:
                    cds_start, cds_end = int(cds_tuple[0]), int(cds_tuple[1])
                    for index, pos in enumerate(targetpos[i]):
                        reg = classify_region(
                            pos, target_end[index], cds_start, cds_end
                        )
                        # There is no 5'UTR when using the annotation file:
                        # strip any "5'UTR" labels; if nothing remains, treat as unclassified.
                        if isinstance(reg, list):
                            reg = [r for r in reg if r != "5'UTR"]
                            reg = reg if len(reg) else None
                        regions.append(reg)
                else:
                    regions = [None] * len(targetpos[i])

        elif input_type == "fasta":
            # look up CDS from Ensembl cds/cdna files by transcript id
            transcript = header
            cds_file = config.cds_file
            cdna_file = config.cdna_file
            df = find_start_end_sites(cds_file, cdna_file, transcript)
            if not df.empty:
                cds_start = int(df["CDS_start"].values[0])
                cds_end = int(df["CDS_end"].values[0])
                for index, pos in enumerate(targetpos[i]):
                    regions.append(
                        classify_region(pos, target_end[index], cds_start, cds_end)
                    )
            else:
                regions = [None] * len(targetpos[i])
        elif input_type == "csv":
            #get the ensembl gene symbol
            print(fullheader)
            gene_symbol = get_gene_symbol(fullheader)
            print(f"gene symbol: {gene_symbol}")
            #find entries in the cds database for this gene symbol
            cds_headers = str(config.cds_file) + '_headers.txt'
            print(cds_headers)
            cds_entries =  find_cds_entries(cds_headers, gene_symbol)
            print(f'cds entries: {cds_entries}')
            #Iterate, check if it's in the variant list, if it's not, bad luck, complain. 
            #define variant list
            variants = find_latest_variants_fasta(tempdir)
            print(variants)
            for cds_variants in cds_entries:
                variant_id = cds_variants.split()[0].lstrip('>')
                if variant_id in variants:
                    print(f'Found valid (some annotated cds) variant: {variant_id}')
                    reference_variant = variant_id
                    break
                else:
                    print(f'No valid (some annotated cds) variant for: {variant_id}')
            #map the target to the reference variant sequence
            reference_variant_sequence = get_cdna(reference_variant)
            results = [align_to_reference_variant(target_sequence, reference_variant_sequence) for target_sequence in targets[0]]
            #it's a tuple of (start in ref, length of match, matched sequence)
            print(results)
            #map the reference variant cds into the cdna
            df = find_start_end_sites(config.cds_file, config.cdna_file, reference_variant)
            print(df)
            if df.empty:
                print('No CDS found for variant')
                regions = [None]*len(probes)
            else:
                cds_start = int(df["CDS_start"].values[0])
                cds_end = int(df["CDS_end"].values[0])
                for probe in results:
                    print(probe[0], probe[0]+probe[1], cds_start, cds_end)
                    regions.append(
                        classify_region(probe[0], probe[0]+probe[1], cds_start, cds_end)
                    )
        print(f"Found Regions: {regions}")
        binding_regions = config.binding_regions
        # Keep only probes whose classified regions intersect the allowed binding regions
        filtered_results = [
            (index, region)
            for index, region in enumerate(regions)
            if region and any(elem in binding_regions for elem in region)
        ]
        if filtered_results:
            filtered_indices, filtered_regions[i] = zip(*filtered_results)
            filtered_regions[i] = list(filtered_regions[i])
        else:
            filtered_indices, filtered_regions[i] = [], []

        # Log region-based drops (before applying the filter)
        prefilter_count = len(regions)
        if prefilter_count:
            allowed_set = set(binding_regions)
            kept_set = set(filtered_indices)
            for idx in range(prefilter_count):
                if idx in kept_set:
                    continue
                reg = regions[idx]
                reason = "unclassified region (no match to CDS/annotation)"
                if reg:
                    if not any(elem in allowed_set for elem in reg):
                        reason = f"region not allowed (region={reg}, allowed={binding_regions})"
                print(
                    f"Dropping (region filter) header={header} idx={idx} reason={reason}; "
                    f"target_pos={targetpos[i][idx] if idx < len(targetpos[i]) else 'NA'}; "
                    f"Tm={(Tm[i][idx] if idx < len(Tm[i]) else 'NA')}; "
                    f"probe_seq={probes[i][idx] if idx < len(probes[i]) else 'NA'}"
                )

        # If nothing passed the region filter for this header, skip to next
        if not filtered_indices:
            continue

        # Apply the region filter to probes, Tm, targetpos, targets
        probes[i] = [probes[i][j] for j in filtered_indices]
        Tm[i] = [Tm[i][j] for j in filtered_indices]
        targetpos[i] = [targetpos[i][j] for j in filtered_indices]
        targets[i] = [targets[i][j] for j in filtered_indices]

        # Group indices by region to select n per region bucket
        region_to_indices = {}
        for index, region in enumerate(filtered_regions[i]):
            region_key = tuple(region) if isinstance(region, list) else region
            if region_key not in region_to_indices:
                region_to_indices[region_key] = []
            region_to_indices[region_key].append(index)

        for region, indices in region_to_indices.items():
            # Guard against stale indices
            indices = [c for c in indices if c < len(probes[i])]
            if not indices:
                continue

            # probes with homopolymers
            wAAAA = [c for c in indices if "AAAA" in probes[i][c]]
            wCCCC = [c for c in indices if "CCCC" in probes[i][c]]
            wGGGG = [c for c in indices if "GGGG" in probes[i][c]]
            wTTTT = [c for c in indices if "TTTT" in probes[i][c]]
            wHomo = set(wAAAA + wCCCC + wGGGG + wTTTT)

            # without homopolymers
            noHomo = list(set(indices) - wHomo)

            # Track per-candidate reasons for potential dropping
            reasons_map = {c: [] for c in indices}
            for c in wHomo:
                reasons_map[c].append(
                    "contains 4-base homopolymer in padlock (AAAA/CCCC/GGGG/TTTT)"
                )

            # enough base complexity (do not consider low complexity at all)
            complexbase = []
            simplebase = []
            for c in indices:
                target = targets[i][c]
                # remove any target has longer than 5 stretch of the same base
                if (
                    "GGGGG" in target
                    or "AAAAA" in target
                    or "CCCCC" in target
                    or "TTTTT" in target
                ):
                    simplebase.append(c)
                    reasons_map[c].append(">=5 identical base run in target")
                else:
                    # Chop the target sequence into substrings of length 10 with a step of 5
                    substring = chopseq(target, 10, 5)

                    # Count the number of unique bases in each substring
                    nbase = [len(set(j)) for j in substring]

                    # If any substring has only 1 or 2 unique bases, it is considered simple
                    if 1 in nbase or 2 in nbase:
                        simplebase.append(c)
                        reasons_map[c].append(
                            "low complexity in 10-mer window (<=2 unique bases)"
                        )
                    else:
                        # Chop the target sequence into substrings of length 2 with a step of 2
                        substring = chopseq(target, 2, 2)

                        # Get the unique substrings
                        unique_substring = list(set(substring))

                        # Count the occurrence of each unique substring
                        ndoublets = [substring.count(i) for i in unique_substring]

                        # If the most common substring appears less than 4 times, it is considered complex
                        if max(ndoublets) < 4:
                            complexbase.append(c)
                        else:
                            # Get the indices of the most common substrings
                            idx = [
                                j
                                for j, tmp in enumerate(ndoublets)
                                if ndoublets[j] == max(ndoublets)
                            ]

                            # Assume the sequence is not simple
                            simple = False

                            # If any of the most common substrings is not a homopolymer and appears 4 times consecutively in the target, it is considered simple
                            for j in idx:
                                if (
                                    unique_substring[j] not in ["AA", "CC", "GG", "TT"]
                                    and unique_substring[j] * 4 in target
                                ):
                                    simple = True
                                    break

                            # If the sequence is simple, add it to the simple base list
                            if simple:
                                simplebase.append(c)
                                reasons_map[c].append(
                                    "repeated non-homopolymer doublet x4 consecutively in target"
                                )
                            else:
                                # Otherwise, add it to the complex base list
                                complexbase.append(c)

            # probes ranking
            primary_targets = list(set(noHomo) & set(complexbase))
            secondary_targets = list(wHomo & set(complexbase))

            # prioritize sequence without homopolymers and no repeated substrings
            random_pruned = set()
            if len(primary_targets) > n:
                pruned = random.sample(primary_targets, len(primary_targets) - n)
                random_pruned = set(pruned)
                deletei = pruned + list(wHomo | set(simplebase))
            elif len(primary_targets) + len(secondary_targets) > n:
                pruned = random.sample(
                    secondary_targets, len(secondary_targets) - n + len(primary_targets)
                )
                random_pruned = set(pruned)
                deletei = pruned + simplebase
            else:
                deletei = (
                    simplebase  # if still not enough, get rid of low-complexity ones
                )

            # Log per-probe drop reasons before deleting
            seen = set()
            deletei.sort(reverse=True)
            for j in deletei:
                if j in seen or j >= len(probes[i]):
                    continue
                seen.add(j)
                reasons = list(reasons_map.get(j, []))
                if j in random_pruned:
                    reasons.append("downselected to n per region")
                if not reasons:
                    reasons = ["unspecified"]
                print(
                    f"Dropping (selection) header={header} region={region} idx={j}; "
                    f"reason={'; '.join(reasons)}; target_pos={targetpos[i][j]}; "
                    f"Tm={Tm[i][j]}; probe_seq={probes[i][j]}"
                )
                del targets[i][j]
                del Tm[i][j]
                del targetpos[i][j]
                del probes[i][j]
                del filtered_regions[i][j]

    return (probes, Tm, targetpos, targets, filtered_regions)


def find_regions(cdna_seq, cds_seq):
    start_cds = cdna_seq.find(cds_seq)
    if start_cds == -1:
        return None  # If the CDS sequence is not found in the cdna sequence
    end_cds = start_cds + len(cds_seq)
    return start_cds, end_cds


# Function to parse a FASTA file and return a dictionary with metadata
def parse_fasta(file_path):
    sequences = {}
    metadata = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
        header = record.description

        # Extract gene_symbol
        gene_symbol = ""
        if "gene_symbol:" in header:
            gene_symbol = header.split("gene_symbol:")[1].split()[0]

        # Extract transcript_biotype
        transcript_biotype = ""
        if "transcript_biotype:" in header:
            transcript_biotype = header.split("transcript_biotype:")[1].split()[0]

        # Extract description
        description = ""
        if "description:" in header:
            description = header.split("description:")[1].strip()

        metadata[record.id] = {
            "gene_symbol": gene_symbol,
            "transcript_biotype": transcript_biotype,
            "description": description,
        }
    return sequences, metadata


def find_start_end_sites(cds_file, cdna_file, transcript):
    """Load ensembl files and find start and end sites of CDS for each transcript.

    Args:
        cds_file (path): Path to the CDS ensemble file.
        cdna_file (path): Path to the cDNA ensembl file.
        transcript (str): Transcript ID to search for.

    Returns:
        DataFrame: DataFrame containing the start and end sites of the CDS for each
        transcript in the list.
    """
    cds_sequences, _ = parse_fasta(cds_file)
    cdna_sequences, cdna_metadata = parse_fasta(cdna_file)

    results = []
    for transcript_id, cdna_seq in cdna_sequences.items():
        if transcript_id in cds_sequences:
            cds_seq = cds_sequences[transcript_id]
            regions = find_regions(cdna_seq, cds_seq)
            if regions:
                start_cds, end_cds = regions
                metadata = cdna_metadata.get(transcript_id, {})
                result = {
                    "Transcript_ID": transcript_id,
                    "5'UTR_start": 0,
                    "CDS_start": start_cds,
                    "CDS_end": end_cds,
                    "3'UTR_end": len(cdna_seq),
                    "gene_symbol": metadata.get("gene_symbol", ""),
                    "transcript_biotype": metadata.get("transcript_biotype", ""),
                    "description": metadata.get("description", ""),
                }
                results.append(result)
    df = pd.DataFrame(results)

    # Filter the dataframe to only include the transcripts in transcript
    df = df[df["Transcript_ID"] == transcript]

    return df


def classify_region(target_start, target_end, cds_start, cds_end, include_5utr=True):
    """Classify probe binding regions based on startpos, endpos, cds_start, and cds_end

    Args:
        target_start (int): Start position of the probe binding region
        target_end (int): End position of the probe binding region
        cds_start (int): Start position of the CDS
        cds_end (int): End position of the CDS
        include_5utr (bool): If False, do not return "5'UTR" (used when annotation file lacks 5'UTR)
    """
    if (
        target_start is None
        or target_end is None
        or cds_start is None
        or cds_end is None
    ):
        return None

    if not include_5utr:
        # Treat pre-CDS-only binding as not classifiable (since there's no 5'UTR in the annotation file)
        if target_end < cds_start:
            return None
        elif (target_start < cds_start) and (target_end <= cds_end):
            return ["CDS"]
        elif (target_start < cds_start) and (target_end > cds_end):
            return ["CDS", "3'UTR"]
        elif (target_start >= cds_start) and (target_end <= cds_end):
            return ["CDS"]
        elif (
            (target_start >= cds_start)
            and (target_start <= cds_end)
            and (target_end > cds_end)
        ):
            return ["CDS", "3'UTR"]
        elif target_start > cds_end:
            return ["3'UTR"]
        else:
            return None

    # include_5utr == True (original behavior)
    if target_end < cds_start:
        return ["5'UTR"]
    elif (target_start < cds_start) and (target_end <= cds_end):
        return ["5'UTR", "CDS"]
    elif (target_start >= cds_start) and (target_end <= cds_end):
        return ["CDS"]
    elif (target_start < cds_start) and (target_end > cds_end):
        return ["5'UTR", "CDS", "3'UTR"]
    elif (
        (target_start >= cds_start)
        and (target_start <= cds_end)
        and (target_end > cds_end)
    ):
        return ["CDS", "3'UTR"]
    elif target_start > cds_end:
        return ["3'UTR"]
    else:
        return None

def align_to_reference_variant(target_sequence, reference_variant_seq):
    # target_sequence = short probe
    # reference_variant_seq = long cDNA
    matcher = SequenceMatcher(None, reference_variant_seq, target_sequence)
    match = matcher.find_longest_match(0, len(reference_variant_seq),
                                       0, len(target_sequence))
    #check that the match fully matches? TODO
    return match.a, match.size, reference_variant_seq[match.a:match.a + match.size]

def get_gene_symbol(header):
    # Split by spaces
    parts = header.split()

    # Find the token that starts with "gene_symbol:"
    gene_symbol_field = next((p for p in parts if p.startswith("gene:")), None)

    # Extract the value after the colon
    if gene_symbol_field:
        gene_symbol = gene_symbol_field.split(":", 1)[1]
    else:
        gene_symbol = None

    return gene_symbol


def find_latest_variants_fasta(tempdir):
    # Look for files like variants_round3.fasta in the given directory
    pattern = os.path.join(tempdir, "*variants_round*.fasta")
    round_files = glob.glob(pattern)

    if round_files:
        print("found some variant files")

        # Extract round numbers using regex and pick the highest one
        def round_num(fpath):
            fname = os.path.basename(fpath)
            m = re.search(r"variants_round(\d+)\.fasta$", fname)
            return int(m.group(1)) if m else -1

        latest = max(round_files, key=round_num)
        print(f"picking latest: {latest}")

        with open(latest, "r") as f:
            lines = f.readlines()
            headers = [line.strip() for line in lines if line.startswith(">")]
            print(headers)
            variants = [h.lstrip(">").split()[0] for h in headers]

        return variants

    # Fallback: look for variants.fasta
    pattern = os.path.join(tempdir, "*variants.fasta")
    plain_files = glob.glob(pattern)  
    try:
        plain_file = plain_files[0]
    except IndexError:
        print("No variant files found")
        return None
    if os.path.exists(plain_file):
        with open(plain_file, "r") as f:
            lines = f.readlines()
            headers = [line.strip() for line in lines if line.startswith(">")]
            print(headers)
            variants = [h.lstrip(">").split()[0] for h in headers]
        return variants

def find_cds_entries(cds_headers_file, gene_symbol):
    entries = []
    with open(cds_headers_file, 'r') as f:
        lines = f.readlines()
        cds_headers = [line.strip() for line in lines]
    for header in cds_headers:
        if f"gene:{gene_symbol}" in header:
            entries.append(header)
    return entries

def get_cdna(transcript_id):
    """
    Return the cDNA sequence (as a string) for a given transcript ID
    from a multi-entry FASTA file.
    """
    for record in SeqIO.parse(config.cdna_file, "fasta"):
        if record.id == transcript_id:
            return str(record.seq)
    raise ValueError(f"Transcript ID {transcript_id} not found in {config.cdna_file}")