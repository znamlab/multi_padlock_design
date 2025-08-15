# retrieve sequence if not given for probe design
# for multi_padlock_design
# Xiaoyan, 2016-8-4


import os
import config
from lib import parmsa
from lib.formatrefseq import fastadb


Acronyms = []
Seq = []
Headers = []


def loaddb(species):
    """Load formatted RefSeq header and sequence files"""
    global Acronyms
    global Headers
    global Seq

    try:
        fastadir = os.path.join(config.BASE_DIR, config.reference_transcriptome)
        if getattr(config, "reference_transcriptome", "").lower() == "refseq":
            # Expect config.fasta_pre_suffix_refseq = (prefix, suffix), config.fasta_filenum_refseq = int
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
                        "Invalid RefSeq config: expected fasta_pre_suffix_refseq=(prefix, suffix) and fasta_filenum_refseq=int"
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

    hits = []
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
    for c, hit in enumerate(hits):
        if len(hit) == 1:  # only one variant
            targets.append(Seq[hit[0]])
            headers.append(Headers[hit[0]])
            basepos.append([0, len(Seq[hit[0]]) - 1])
            variants.append(Headers[hit[0]][1:].split(".", 1)[0])
        else:
            msa.append(genes[c])
            tempheader = []
            file = dirname + "/" + genes[c] + "_variants.fasta"
            with open(file, "w") as f:
                for multi in hit:
                    f.write("%s\n" % Headers[multi])
                    f.write("%s\n\n" % Seq[multi])
                    tempheader.append(Headers[multi])
            headersMSA.append(tempheader)
            variants.append([i[1:].split(".", 1)[0] for i in tempheader])

    # run multiple sequence alignment if more than one variants are found for one gene
    if len(msa):
        tempout = parmsa.continuemsa(dirname, msa)
        print("MSA finished.")
        for c, name in enumerate(tempout[0]):
            if len(tempout[1][c]):
                tempheader = headersMSA[c]
                temp = [
                    i for i in tempheader if name in i
                ]  # first sequence in ClustalW output
                headers.append(temp[0])
                basepos.append(tempout[1][c])
                targets.append(tempout[2][c])
            else:
                nocommon.append(msa[c])

    return headers, basepos, targets, msa, nocommon, variants
