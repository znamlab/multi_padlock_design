#!/bin/bash
#SBATCH --job-name=blast_${1%.*}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --partition=ncpu
#SBATCH -e /nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/probe_checking/probe_${1%.*}.err

#tmux
sleep 3

ml Anaconda3/2022.05
ml BLAST+

source /camp/apps/eb/software/Anaconda/conda.env.sh
conda activate iss-preprocess
cd /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/padlock_checking

REFSEQ_DIR=/nemo/lab/znamenskiyp/home/shared/resources/refseq
SCRATCH_DIR=/nemo/lab/znamenskiyp/scratch

#blastn -query $INPUT -db "${SCRATCH_DIR}/test_refseq/mouse.transcriptome" -outfmt 10 -out "${INPUT%.*}_blast.out" -word_size 7 -strand plus
blastn -query $INPUT -db "${SCRATCH_DIR}/test_refseq/mouse.transcriptome" -outfmt "10 std qseq sseq" -out "${INPUT%.*}_blast.out" -word_size 7 -strand plus

