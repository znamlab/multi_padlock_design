import os
from pathlib import Path

BASE_DIR = os.path.abspath("/nemo/lab/znamenskiyp/home/shared/resources/")

# Directory to save split input files into
# e.g., an input file "Olfrs.fa" will be split into:
#   /nemo/lab/znamenskiyp/scratch/Olfrs/<gene>.fasta
split_input = Path("/nemo/lab/znamenskiyp/scratch")

# maximum number of processes for parallel ClustalW
max_MSA_processes = 32
# maximum number of processes for parallel Blast
max_Blast_processes = 32

# Database to use for BLAST specificity checking
reference_transcriptome = "ensembl"
species = "mouse"

# Refseq database naming convention
fasta_filenum_refseq = 3  # Change this for different species, 3 for mouse
fasta_pre_suffix_refseq = ("mouse.", ".rna.fna")

# Ensembl database filenames
# cds_file = r"Mus_musculus.GRCm39.cds.all.fa"
# cdna_file = r"Mus_musculus.GRCm39.cdna.all.fa"
cds_file = os.path.join(
    BASE_DIR, reference_transcriptome, r"Mus_musculus.GRCm39.cds.all.fa"
)
cdna_file = os.path.join(
    BASE_DIR, reference_transcriptome, r"Mus_musculus.GRCm39.cdna.all.fa"
)
gencode_file = os.path.join(BASE_DIR, "gencode", r"gencode.vM38.basic.annotation.gtf")
# extra_files = [r"olfr_transcripts.fa"]

# Aggregated BLAST file
blast_db_file = os.path.join(
    BASE_DIR, reference_transcriptome, f"{species}.selected.fasta"
)

# mRNA regions to target for probe binding
binding_regions = ["CDS", "3'UTR"]
# Use melting temperature for BLAST specificity checking
specificity_by_tm = True

# OLFR annotation file for custom olfr probe design
# annotation_file = os.path.join(
#    BASE_DIR, "olfr_annotations", "olfr_consensus_cds_3utr_annotations.csv"
# )
