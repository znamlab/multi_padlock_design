import os
from pathlib import Path

BASE_DIR = os.path.abspath("/nemo/lab/znamenskiyp/home/shared/resources/")
# Resolve the repository root by walking up from this config module.

REPO_ROOT = Path(__file__).resolve().parent.parent
# Directory to save split input files into
# e.g., an input file "Olfrs.fa" will be split into:
#   /nemo/lab/znamenskiyp/scratch/Olfrs/<gene>.fasta
# split_input = Path("/nemo/lab/znamenskiyp/scratch")
split_input = REPO_ROOT / "output"

# Directory that stores per-job probe design outputs
# probe_output_root = Path("/nemo/lab/znamenskiyp/scratch")
probe_output_root = REPO_ROOT / "output"

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

# Per-probe (padlock-specific) identifiers
# When enabled, each padlock within a gene carries a short identifier inserted
# between the barcode and the downstream linker, so individual probes of the same
# gene can be told apart in the readout. Applied to the final random-subset
# deliverable only (requires a per-gene probe count, i.e. "Number of probes per
# gene" must be set).
use_padlock_identifiers = False
padlock_identifier_file = os.path.join(
    REPO_ROOT, "multi_padlock_design", "barcodes", "padlock_identifiers_3bp.txt"
)
# Minimum Hamming distance defining the error-robust tier of the identifier set.
# The leading codes meeting this distance are used first; codes beyond the robust
# tier are still unique but only distance-1 (a warning is emitted when used).
padlock_identifier_min_distance = 2

# Default CLI parameters for slurm probe design submissions
probe_design_defaults = {
    "species": species,
    "arm_len": 20,
    "total_len": 40,
    "interval": 5,
    "tm_low": 60,
    "tm_high": 78,
    "n_probes": 100,
}

# OLFR annotation file for custom olfr probe design
# annotation_file = os.path.join(
#    BASE_DIR, "olfr_annotations", "olfr_consensus_cds_3utr_annotations.csv"
# )
