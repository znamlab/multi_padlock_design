# for multi_padlock_design
# Xiaoyan, 2018

import random
from lib.screenseq import chopseq
from Bio import SeqIO
import config

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


def selectprobes(n, finals, headers, armlength):
    """Prioritize probes with no homopolymer sequences and choose randomly n candidates per region

    Args:
        n (int): number of probes per region
        finals (list): list of final probes
        headers (list): list of headers
        armlength (int): arm length
    
    Returns:
        probes (list): list of final probes
        Tm (list): list of Tm values
        targetpos (list): list of target positions
        targets (list): list of targets
        filtered_regions (list): list of binding regions
    """
    print(n)
    print(finals)
    print(headers)
    print(armlength)
    probes = finals[0]
    Tm = finals[1]
    targetpos = finals[2]
    targets = finals[3]
    filtered_regions = [0] * len(headers)
    for i, header in enumerate(headers):
        # Find transcript region of probe binding site
        transcript = header[1:].split(" ")[0]
        cds_file = config.cds_file
        cdna_file = config.cdna_file
        df = find_start_end_sites(cds_file, cdna_file, transcript)
        target_end = [target + (1 + (2 * armlength)) for target in targetpos[i]]
        regions = []
        for index, pos in enumerate(targetpos[i]):
            regions.append(classify_region(pos, target_end[index], df["CDS_start"].values, df["CDS_end"].values))
        binding_regions = config.binding_regions
        # If the probe binding site is not in the list of binding_regions, remove the probe
        filtered_results = [(index, region) for index, region in enumerate(regions) if any(elem in binding_regions for elem in region)]
        if filtered_results:
            filtered_indices, filtered_regions[i] = zip(*filtered_results)
            filtered_regions[i] = list(filtered_regions[i])
        else:
            filtered_indices, filtered_regions[i] = [], []

        # Check if filtered_indices is empty
        if not filtered_indices:
            continue  # Skip this iteration if there are no valid indices

        # Filter probes, Tm, targetpos, targets based on filtered indices
        probes[i] = [probes[i][j] for j in filtered_indices]
        Tm[i] = [Tm[i][j] for j in filtered_indices]
        targetpos[i] = [targetpos[i][j] for j in filtered_indices]
        targets[i] = [targets[i][j] for j in filtered_indices]

        # Group probes by region
        region_to_indices = {}
        for index, region in enumerate(filtered_regions[i]):
            region_key = tuple(region) if isinstance(region, list) else region
            if region_key not in region_to_indices:
                region_to_indices[region_key] = []
            region_to_indices[region_key].append(index)

        for region, indices in region_to_indices.items():
            # Check if indices are within the range of probes[i]
            indices = [c for c in indices if c < len(probes[i])]
            if not indices:
                continue  # Skip if no valid indices

            # probes with homopolymers
            wAAAA = [c for c in indices if "AAAA" in probes[i][c]]
            wCCCC = [c for c in indices if "CCCC" in probes[i][c]]
            wGGGG = [c for c in indices if "GGGG" in probes[i][c]]
            wTTTT = [c for c in indices if "TTTT" in probes[i][c]]
            wHomo = set(wAAAA + wCCCC + wGGGG + wTTTT)

            # without homopolymers
            noHomo = list(set(indices) - wHomo)

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
                else:
                    # Chop the target sequence into substrings of length 10 with a step of 5
                    substring = chopseq(target, 10, 5)

                    # Count the number of unique bases in each substring
                    nbase = [len(set(j)) for j in substring]

                    # If any substring has only 1 or 2 unique bases, it is considered simple
                    if 1 in nbase or 2 in nbase:
                        simplebase.append(c)
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
                            else:
                                # Otherwise, add it to the complex base list
                                complexbase.append(c)

            # probes ranking
            primary_targets = list(set(noHomo) & set(complexbase))
            secondary_targets = list(wHomo & set(complexbase))

            # prioritize sequence without homopolymers and no repeated substrings
            if len(primary_targets) > n:
                deletei = random.sample(primary_targets, len(primary_targets) - n)
                deletei = deletei + list(wHomo | set(simplebase))
            elif len(primary_targets) + len(secondary_targets) > n:
                deletei = random.sample(
                    secondary_targets, len(secondary_targets) - n + len(primary_targets)
                )
                deletei = deletei + simplebase
            else:
                deletei = simplebase  # if still not enough, get rid of low-complexity ones

            deletei.sort(reverse=True)
            for j in deletei:
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
            "description": description
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
    import pandas as pd
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
                    "description": metadata.get("description", "")
                }
                results.append(result)
    df = pd.DataFrame(results)

    # Filter the dataframe to only include the transcripts in transcript
    df = df[df["Transcript_ID"] == transcript]
    
    return df

def classify_region(target_start, target_end, cds_start, cds_end):
    """Classify probe binding regions based on startpos, endpos, cds_start, and cds_end

    Args:
        target_start (int): Start position of the probe binding region
        target_end (int): End position of the probe binding region
        cds_start (int): Start position of the CDS
        cds_end (int): End position of the CDS

    Returns:
        list: List of regions that the probe binding region falls into
    """

    if target_start is None or target_end is None:
        return None

    if target_end < cds_start:
        return ["5'UTR"]
    elif (target_start < cds_start) and (target_end <= cds_end):
        return ["5'UTR", "CDS"]
    elif (target_start >= cds_start) and (target_end <= cds_end):
        return ["CDS"]
    elif (target_start < cds_start) and (target_end > cds_end):
        return ["5'UTR", "CDS", "3'UTR"]
    elif (target_start >= cds_start) and (target_start <= cds_end) and (target_end > cds_end):
        return ["CDS", "3'UTR"]
    elif target_start > cds_end:
        return ["3'UTR"]
    else:
        return None