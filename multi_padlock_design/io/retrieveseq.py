# retrieve sequence if not given for probe design
# for multi_padlock_design
# Xiaoyan, 2016-8-4


import os
import re

import pandas as pd
from Bio import SeqIO

import multi_padlock_design.config as config
from multi_padlock_design.io.formatrefseq import fastadb
from multi_padlock_design.io.gene_utils import get_gene_synonyms
from multi_padlock_design.msa import parmsa

Acronyms = []
Seq = []
Headers = []
Gencode = []


def loaddb(species):
    """Load formatted RefSeq header and sequence files"""
    global Acronyms
    global Headers
    global Seq
    global Gencode

    try:
        fastadir = os.path.join(config.BASE_DIR, config.reference_transcriptome)
        if getattr(config, "reference_transcriptome", "").lower() == "refseq":
            # Expect config.fasta_pre_suffix_refseq = (prefix, suffix),
            # config.fasta_filenum_refseq = int
            pre_suf = getattr(config, "fasta_pre_suffix_refseq", None)
            filenum = getattr(config, "fasta_filenum_refseq", None)
            acr_path = os.path.join(fastadir, species + ".acronymheaders.txt")
            if not os.path.isfile(acr_path):
                print("Processing fasta database files..")
                if (
                    isinstance(pre_suf, (list, tuple))
                    and len(pre_suf) == 2
                    and isinstance(filenum, int)
                ):
                    fastadb(
                        fastadir,
                        (pre_suf[0], pre_suf[1], filenum),
                        species,
                        source="refseq",
                        keep_nm_nr_only=True,
                    )
                else:
                    raise ValueError(
                        "Invalid RefSeq config: expected fasta_pre_suffix_refseq="
                        "(prefix, suffix) and fasta_filenum_refseq=int"
                    )

        elif getattr(config, "reference_transcriptome", "").lower() == "ensembl":
            files = [config.cdna_file]
            if hasattr(config, "extra_files") and config.extra_files:
                files += config.extra_files
            acr_path = os.path.join(fastadir, species + ".acronymheaders.txt")
            if not os.path.isfile(acr_path):
                print("Processing fasta database files..")
                fastadb(
                    fastadir,
                    files,
                    species,
                    source="ensembl",
                    keep_nm_nr_only=False,
                )

        # load database files
        with open(os.path.join(fastadir, species + ".acronymheaders.txt"), "r") as f:
            Acronyms = [line.rstrip("\n") for line in f]
        with open(os.path.join(fastadir, species + ".selectedheaders.txt"), "r") as f:
            Headers = [line.rstrip("\n") for line in f]
        with open(os.path.join(fastadir, species + ".selectedseqs.txt"), "r") as f:
            Seq = [line.rstrip("\n") for line in f]
        Gencode = load_gencode()

    except FileNotFoundError:
        print("Cannot load fasta database due to mismatch in species.")


def querygenes(genes, species):
    """Find genes from a database"""
    # species check performed in func. keyboardinput
    global Acronyms
    global Headers
    global Seq

    # load specified database
    loaddb(species)

    def get_ensembl_id(header):
        match = re.search(r"gene:([^ ]+)", header)
        if match:
            return match.group(1)
        return None

    hits = []
    if config.reference_transcriptome == "ensembl":
        print("filtering by GENCODE Basic annotation")
        for gene in genes:
            if gene not in Acronyms:
                translated_id = get_gene_synonyms(
                    config.species, gene
                )  # Correct gene name if possible
                if translated_id not in Acronyms:
                    hits.append([])
                    continue
                else:
                    gene = translated_id
            # get ensembl ID
            ensembl_ID = [
                get_ensembl_id(header) for header in Headers if gene in header
            ][0]
            print(f"Querying for {ensembl_ID}")
            # filter by Gencode basic annotation
            allowed_variants = Gencode[Gencode["gene_id"] == ensembl_ID][
                "transcript_id"
            ].tolist()
            supporting_evidence = Gencode[Gencode["gene_id"] == ensembl_ID][
                "transcript_support_level"
            ].tolist()
            print(
                f"allowed variants: {allowed_variants},"
                f" supporting evidence: {supporting_evidence} (1 is best)"
            )
            # now form the list of all variants and provide a mask
            hit = [
                c for c, header in enumerate(Headers) if ensembl_ID in header
            ]  # variant search now happens via ensembl ID
            # filter hits by allowed variants
            hit_masked = [(Headers[h].split()[0][1:] in allowed_variants) for h in hit]
            # filter hits by allowed
            hits.append(hit)

        return (hits, hit_masked)

    else:
        for gene in genes:
            if gene not in Acronyms:
                hits.append([])
            else:
                hit = [c for c, header in enumerate(Acronyms) if header == gene]
                hits.append(hit)

    return hits


def findseq(genes, hits, dirname):
    """Find target sequences, if multiple variants are found, call MSA"""
    global Acronyms
    global Headers
    global Seq

    targets = []
    headers = []
    basepos = []
    msa = []
    nocommon = []
    variants = []

    headersMSA = []

    if config.reference_transcriptome == "ensembl":
        hits, hitmask = hits[0], hits[1]

    for c, hit in enumerate(hits):
        if len(hit) == 1:  # only one variant
            targets.append(Seq[hit[0]])
            headers.append(Headers[hit[0]])
            basepos.append([0, len(Seq[hit[0]]) - 1])
            variants.append(Headers[hit[0]][1:].split(".", 1)[0])
            file = (
                dirname + "/" + genes[c] + "_variants.fasta"
            )  # makes it easier to map properly the target
            # if all sets of variants are written down
            with open(file, "w") as f:
                f.write("%s\n" % headers[0])
                f.write("%s\n\n" % targets[0])
            full_variants = (
                variants  # no filtering happened, since there's only one variant
            )
        else:  # more than ona variant
            msa.append(genes[c])
            tempheader = []
            file = dirname + "/" + genes[c] + "_variants.fasta"
            with open(file, "w") as f:
                full_variants = []
                for multi in hit:
                    f.write("%s\n" % Headers[multi])
                    f.write("%s\n\n" % Seq[multi])
                    full_variants.append(Headers[multi][1:].split(".", 1)[0])
            if (
                config.reference_transcriptome == "ensembl"
            ):  # only keep variants that pass the GENCODE mask
                hit = [h for i, h in enumerate(hit) if hitmask[i]]
                if len(hit) == 1:
                    msa = []  # we CLUSTAL as a function of the filtered variants
                    targets.append(Seq[hit[0]])
                    headers.append(Headers[hit[0]])
                    basepos.append([0, len(Seq[hit[0]]) - 1])
                file = dirname + "/" + genes[c] + "_allowed_variants.fasta"
                with open(file, "w") as f:
                    for multi in hit:
                        f.write("%s\n" % Headers[multi])
                        f.write("%s\n\n" % Seq[multi])
                        tempheader.append(Headers[multi])
                headersMSA.append(tempheader)
                variants.append([i[1:].split(".", 1)[0] for i in tempheader])
            else:  # Clustal on all variants, no filter
                for multi in hit:
                    tempheader.append(Headers[multi])
                headersMSA.append(tempheader)
                variants.append([i[1:].split(".", 1)[0] for i in tempheader])

    # run multiple sequence alignment if more than one variant is found for any gene
    if len(msa):
        # process each gene independently
        for gi, gene in enumerate(msa):
            # headers for all variants of this gene (written earlier)
            tempheader = headersMSA[gi]

            round_no = 1  # first run uses *_variants.fasta / .aln / .dnd

            def fasta_for_round(r: int) -> str:
                if config.reference_transcriptome == "ensembl":
                    return (
                        f"{dirname}/{gene}_allowed_variants.fasta"
                        if r == 1
                        else f"{dirname}/{gene}_allowed_variants_round{r}.fasta"
                    )
                return (
                    f"{dirname}/{gene}_variants.fasta"
                    if r == 1
                    else f"{dirname}/{gene}_variants_round{r}.fasta"
                )

            def dnd_for_round(r: int) -> str:
                if config.reference_transcriptome == "ensembl":
                    return (
                        f"{dirname}/{gene}_allowed_variants.dnd"
                        if r == 1
                        else f"{dirname}/{gene}_allowed_variants_round{r}.dnd"
                    )
                return (
                    f"{dirname}/{gene}_variants.dnd"
                    if r == 1
                    else f"{dirname}/{gene}_variants_round{r}.dnd"
                )

            # how many variants do we currently have?
            n_left = sum(1 for _ in SeqIO.parse(fasta_for_round(round_no), "fasta"))

            while True:
                # run MSA for this gene only
                # (pass [gene] so continuemsa handles one job)
                out = parmsa.continuemsa(
                    dirname,
                    [gene],
                    round=None if round_no == 1 else round_no,
                    reset=True,
                )
                # out = (Names, BasePos, Seqs)
                # each is length 1 because we passed [gene]
                name_c, basepos_c, seqs_c = out[0][0], out[1][0], out[2][0]

                if len(basepos_c):  # consensus found
                    # pick the header that matches the first sequence used by ClustalW
                    # (the code relies on matching by substring of the chosen name)
                    chosen = [h for h in tempheader if name_c in h]
                    if not chosen:
                        # fallback: just take the first header if matching fails
                        chosen = [tempheader[0]]
                    headers.append(chosen[0])
                    basepos.append(basepos_c)
                    targets.append(seqs_c)
                    break

                # no consensus: stop if ≤ 2 variants remain
                if n_left <= 2:
                    nocommon.append(gene)
                    break

                # drop one outgroup and try again
                print(
                    f"[{gene}] No consensus, removing outgroup"
                    f" (round {round_no} -> {round_no+1})"
                )
                treefile = dnd_for_round(round_no)
                outgroup = parmsa.find_outgroup(treefile)

                in_fa = fasta_for_round(round_no)
                round_no += 1
                parmsa.remove_outgroup(in_fa, outgroup, round=round_no)

                # update remaining variants
                n_left = sum(1 for _ in SeqIO.parse(fasta_for_round(round_no), "fasta"))

    if config.reference_transcriptome == "ensembl":
        # update variants to blast against all variants later
        variants = full_variants

    print("MSA finished across genes.")

    return headers, basepos, targets, msa, nocommon, variants


def load_gencode():
    def extract_features(line, feature):
        match = re.search(rf'{feature} "([^"]+)"', line)
        if match:
            return match.group(1)
        return None

    gene_ids = []
    transcript_ids = []
    transcript_support_levels = []

    with open(config.gencode_file, "r") as f:
        for line in f:
            if line.startswith("#") or line.split("\t")[2] != "transcript":
                continue  # skip header lines and gene lines
            line = line.strip()
            gene_id = extract_features(line, "gene_id")
            transcript_id = extract_features(line, "transcript_id")
            transcript_support_level = extract_features(
                line, "transcript_support_level"
            )
            if gene_id or transcript_id:
                gene_ids.append(gene_id)
                transcript_ids.append(transcript_id)
                transcript_support_levels.append(transcript_support_level)

    gencode = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "transcript_id": transcript_ids,
            "transcript_support_level": transcript_support_levels,
        }
    )
    return gencode
