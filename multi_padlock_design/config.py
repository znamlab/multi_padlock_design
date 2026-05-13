import os
from pathlib import Path

BASE_DIR = os.path.abspath("/nemo/lab/znamenskiyp/home/shared/resources/")
# Resolve the repository root by walking up from this config module.
REPO_ROOT = Path(__file__).resolve().parent.parent
# Directory to save split input files into
# e.g., an input file "Olfrs.fa" will be split into:
#   /nemo/lab/znamenskiyp/scratch/Olfrs/<gene>.fasta
split_input = Path("/nemo/lab/znamenskiyp/scratch")

# Directory that stores per-job probe design outputs
probe_output_root = Path("/nemo/lab/znamenskiyp/scratch")

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
    BASE_DIR, reference_transcriptome, r"Mus_musculus.GRCm39.release-115.cds.all.fa"
)
cdna_file = os.path.join(
    BASE_DIR, reference_transcriptome, r"Mus_musculus.GRCm39.release-115.cdna.all.fa"
)
gencode_file = os.path.join(BASE_DIR, "gencode", r"gencode.vM38.basic.annotation.gtf")

# cds_file = os.path.join(
#     BASE_DIR, reference_transcriptome, r"Homo_sapiens.GRCh38.release-115.cds.all.fa"
# )
# cdna_file = os.path.join(
#     BASE_DIR, reference_transcriptome, r"Homo_sapiens.GRCh38.release-115.cdna.all.fa"
# )
# gencode_file = os.path.join(BASE_DIR, "gencode", r"gencode.v49.basic.annotation.gtf")
# extra_files = [r"olfr_transcripts.fa"]

# Aggregated BLAST file
blast_db_file = os.path.join(
    BASE_DIR, reference_transcriptome, f"{species}.selected.fasta"
)

# mRNA regions to target for probe binding
binding_regions = ["CDS", "3'UTR"]
# Use melting temperature for BLAST specificity checking
specificity_by_tm = True

# Default CLI parameters for slurm probe design submissions
probe_design_defaults = {
    "species": species,
    "arm_len": 20,
    "total_len": 40,
    "interval": 5,
    "tm_low": 60,
    "tm_high": 78,
    "n_probes": 20,
}

# OLFR annotation file for custom olfr probe design
# annotation_file = os.path.join(
#    BASE_DIR, "olfr_annotations", "olfr_consensus_cds_3utr_annotations.csv"
# )
