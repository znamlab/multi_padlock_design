import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# maximum number of processes for parallel ClustalW
max_MSA_processes = 32
# maximum number of processes for parallel Blast
max_Blast_processes = 32

reference_transcriptome = "ensembl"
species = "mouse"

# mouse refseq database
fasta_filenum_refseq = 3
fasta_pre_suffix_refseq = ("mouse.", ".rna.fna")

# mouse ensemble database
cds_file = r"Mus_musculus.GRCm39.cds.all.fa"
cdna_file = r"Mus_musculus.GRCm39.cdna.all.fa"
# extra_files = [r"olfr_transcripts.fa"]

# Aggregated BLAST file
blast_db_file = os.path.join(
    BASE_DIR, reference_transcriptome, f"{species}.selected.fasta"
)


annotation_file = os.path.join(
    BASE_DIR, "data", "olfr_consensus_cds_3utr_annotations.csv"
)
binding_regions = ["CDS", "3'UTR"]
specificity_by_tm = True
