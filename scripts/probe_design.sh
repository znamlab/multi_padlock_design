#!/bin/bash
#SBATCH --job-name=probe_design
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --partition=ncpu
# Command-line --output from the submitter will override this if provided:
#SBATCH --output=/nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/%x/%x_%j.out

set -euo pipefail

# Expected from sbatch --export:
#   INPUT       absolute path (/path/to/foo.fasta or /path/to/foo.csv)
#   PARENT      basename of input directory
#   INPUT_TYPE  'fasta' or 'csv'

CURRENT_NODE=$(hostname)

resubmit_job() {
  echo "Resubmitting job excluding node: $CURRENT_NODE"
  sbatch --export=ALL --exclude="$CURRENT_NODE" "$0"
  exit 0
}

if command -v ml &>/dev/null; then
  ml purge
  ml Anaconda3/2024.10-1
  ml Clustal-Omega
  ml BLAST+
  ml ClustalW2
else
  echo "Module system is not available on node: $CURRENT_NODE"
  resubmit_job
fi

eval "$(/camp/apps/eb/software/Anaconda3/2024.10-1/bin/conda shell.bash hook)"
conda activate iss-preprocess

cd /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design

FILE_PATH="$INPUT"                              # already absolute
FILE_BASE="$(basename "$FILE_PATH")"
FILE_STEM="${FILE_BASE%.*}"                     # drop .csv/.fasta

OUT_DIR="/nemo/lab/znamenskiyp/scratch/${PARENT}/${FILE_STEM}"
mkdir -p "${OUT_DIR}"

echo "Input file: ${FILE_PATH}"
echo "Mode: ${INPUT_TYPE}"
echo "Output dir: ${OUT_DIR}"
which python

# Common parameters
SPECIES="mouse"
ARM_LEN="20"
INTERVAL="5"
TM_LOW="60"
TM_HIGH="78"
N_PROBES="20"

# Feed probedesign.py prompts
if [[ "${INPUT_TYPE,,}" == "csv" ]]; then
  python probedesign.py <<- ENDOF
${OUT_DIR}
${SPECIES}
${FILE_PATH}
${ARM_LEN}
${INTERVAL}
${TM_LOW}
${TM_HIGH}
${N_PROBES}
ENDOF
elif [[ "${INPUT_TYPE,,}" == "fasta" ]]; then
  python probedesign.py <<- ENDOF
${OUT_DIR}
${SPECIES}

${FILE_PATH}
${ARM_LEN}
${INTERVAL}
${TM_LOW}
${TM_HIGH}
${N_PROBES}
ENDOF
else
  echo "Unknown INPUT_TYPE: ${INPUT_TYPE}" >&2
  exit 1
fi
