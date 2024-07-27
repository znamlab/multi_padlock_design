#!/bin/bash
#SBATCH --job-name=makeblastdb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --partition=ncpu
#SBATCH --mail-type=ALL
#SBATCH -o /nemo/lab/znamenskiyp/home/users/cypranc/slurm_logs/format_refseq_%j.out
#SBATCH -e /nemo/lab/znamenskiyp/home/users/cypranc/slurm_logs/format_refseq_%j.err
#SBATCH --mail-user=cypranc@crick.ac.uk


#tmux
sleep 3

ml Anaconda3/2022.05
ml Clustal-Omega
ml BLAST+
ml ClustalW2

source /camp/apps/eb/software/Anaconda/conda.env.sh
conda activate base
cd /nemo/lab/znamenskiyp/home/users/cypranc/multi_padlock_design

REFSEQ_DIR=/nemo/lab/znamenskiyp/home/shared/resources/refseq
SCRATCH_DIR=/nemo/lab/znamenskiyp/scratch/test_refseq

for i in  $(seq 3)
do
	echo "Making BLAST DB for mouse${i}"
	makeblastdb -in ${REFSEQ_DIR}/mouse.${i}.rna.fna -dbtype nucl -out /nemo/lab/znamenskiyp/scratch/test_refseq
done

blastdb_aliastool -dblist "${REFSEQ_DIR}/mouse.1.rna.fna ${REFSEQ_DIR}/mouse.2.rna.fna ${REFSEQ_DIR}/mouse.3.rna.fna" -dbtype nucl -out ${SCRATCH_DIR}/mouse.transcriptome -title "mouse_transcriptome"
