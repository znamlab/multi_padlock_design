#!/bin/bash
#SBATCH --job-name=par_blast_queries
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --partition=ncpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=becalia@crick.ac.uk

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: sbatch parallel_blast_list.sh <fasta_list.txt> [--skip-existing] [LOG_DIR=/path/to/logs]"
  exit 1
fi

FASTA_LIST="$1"; shift

SKIP_EXISTING=0
LOG_DIR="/nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/probe_checking"

for arg in "$@"; do
  case "$arg" in
    --skip-existing) SKIP_EXISTING=1 ;;
    LOG_DIR=*) LOG_DIR="${arg#LOG_DIR=}" ;;
    *) echo "Unknown argument: $arg" >&2; exit 1 ;;
  esac
done

if [ ! -s "$FASTA_LIST" ]; then
  echo "FASTA list not found or empty: $FASTA_LIST"
  exit 1
fi

mkdir -p "$LOG_DIR"

module purge
module load Anaconda3
module load BLAST+

source /camp/apps/eb/software/Anaconda/conda.env.sh
conda activate iss-preprocess || { echo "Failed to activate conda environment"; exit 1; }

echo "Submitting BLAST jobs for entries in: $FASTA_LIST"
echo "Skip existing outputs: $SKIP_EXISTING"
echo "Logs: $LOG_DIR"

while IFS= read -r fasta; do
  [ -z "$fasta" ] && continue
  if [ ! -f "$fasta" ]; then
    echo "Missing FASTA file (skipping): $fasta"
    continue
  fi
  base=$(basename "$fasta")
  stem="${base%.*}"
  out_candidate="${fasta%_query.fasta}_query_blast.out"

  if [ $SKIP_EXISTING -eq 1 ] && [ -s "$out_candidate" ]; then
    echo "Exists, skipping: $out_candidate"
    continue
  fi

  echo "Launching job for $fasta"
  sbatch \
    --job-name="blast_${stem}" \
    -e "${LOG_DIR}/probe_${stem}.err" \
    --output=/dev/null \
    --export=INPUT="$fasta" \
    blast_query.sh
done < "$FASTA_LIST"

echo "Submission loop complete."