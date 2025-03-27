# maximum number of processes for parallel ClustalW
max_MSA_processes = 32

# maximum number of processes for parallel Blast
max_Blast_processes = 32

# mouse database
fastadir_mouse = r"/nemo/lab/znamenskiyp/home/shared/resources/refseq"
fasta_filenum_mouse = 2
fasta_pre_suffix_mouse = ("mouse.", ".rna.fna")

# human database
fastadir_human = (
    r"/nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/refseq"
)
fasta_filenum_human = 7
fasta_pre_suffix_human = ("human.", ".rna.fna")

cds_file = r"/nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/refseq/Mus_musculus.GRCm39.cds.all.fa"
cdna_file = r"/nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/refseq/Mus_musculus.GRCm39.cdna.all.fa"

binding_regions = ["5'UTR", "CDS", "3'UTR"]
specificity_by_tm = True
