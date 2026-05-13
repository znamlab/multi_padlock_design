#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --partition=ncpu

#tmux
sleep 3

. ~/.bashrc
ml purge
ml Anaconda3
ml BLAST+

source /camp/apps/eb/software/Anaconda/conda.env.sh
conda activate iss-preprocess
cd /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/padlock_checking

# Build transcriptome DB path using config (mirrors logic in parblast.newblast)
DB_PATH=$(python - <<'PY'
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))  # ensure config import
import config
print(os.path.join(config.BASE_DIR, config.reference_transcriptome, config.species + ".transcriptome"))
PY
)

SCRATCH_DIR=/nemo/lab/znamenskiyp/scratch

blastn -query $INPUT -db "$DB_PATH" -outfmt "10 std qseq sseq" -out "${INPUT%.*}_blast.out" -word_size 7 -strand plus
