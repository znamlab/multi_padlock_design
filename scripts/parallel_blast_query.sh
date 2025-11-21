#!/bin/bash
#SBATCH --job-name=par_blast_queries
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
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

# Added: argument parsing & FASTA generation
PADLOCK_CSV="$1"
OUTPUT_DIR="${2:-/nemo/lab/znamenskiyp/scratch/padlock_queries}"
mkdir -p "$OUTPUT_DIR"

if [ -z "$PADLOCK_CSV" ]; then
  echo "Usage: sbatch parallel_blast_query.sh <padlocks.csv> [output_dir]"
  exit 1
fi
if [ ! -f "$PADLOCK_CSV" ]; then
  echo "Input CSV not found: $PADLOCK_CSV"
  exit 1
fi

echo "Generating FASTA files from ${PADLOCK_CSV} into ${OUTPUT_DIR}"
python padlocks_to_fasta.py --input "$PADLOCK_CSV" --output "$OUTPUT_DIR"

FASTA_LIST="${OUTPUT_DIR}/fasta_list.txt"
if [ ! -s "$FASTA_LIST" ]; then
  echo "FASTA list not found or empty: $FASTA_LIST"
  exit 1
fi

LOG_DIR="/nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/probe_checking"
mkdir -p "$LOG_DIR"

# Launch jobs for each FASTA file in the list
cat "$FASTA_LIST" | while read -r line
do
   [ -z "$line" ] && continue
   base=$(basename "$line")
   stem="${base%.*}"
   echo "Launching job for ${line}"
   sbatch \
     --job-name="blast_${stem}" \
     -e "${LOG_DIR}/probe_${stem}.err" \
     --output=/dev/null \
     --export=INPUT="${line}" \
     blast_query.sh || { echo "Failed to submit job for ${line}"; exit 1; }
done
