#!/bin/bash
#SBATCH --job-name=par_blast_queries
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --partition=ncpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=becalia@crick.ac.uk

# Exit on any error
set -e

# Load modules
module purge
module load Anaconda3
module load BLAST+

# Activate conda environment
source /camp/apps/eb/software/Anaconda/conda.env.sh
conda activate iss-preprocess || { echo "Failed to activate conda environment"; exit 1; }

# Navigate to working directory
cd /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/padlock_checking || { echo "Failed to change directory"; exit 1; }

# Set directory variables
REFSEQ_DIR=/nemo/lab/znamenskiyp/home/shared/resources/refseq
FASTA_LIST=/nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/padlock_checking/fasta_list.txt

# Check if FASTA_LIST is set and not empty
if [ -z "$FASTA_LIST" ]; then
  echo "FASTA_LIST variable is not set or empty"
  exit 1
fi

# Launch jobs for each FASTA file in the list
cat "$FASTA_LIST" | while read -r line
do
   echo "Launching job for ${line}"
   sbatch --export=INPUT="${line}" blast_query.sh || { echo "Failed to submit job for ${line}"; exit 1; }
done
