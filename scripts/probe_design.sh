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
#   SPECIES     (optional) species identifier
#   ARM_LEN     (optional) probe arm length
#   TOTAL_LEN   (optional) full probe length
#   INTERVAL    (optional) spacing interval
#   TM_LOW      (optional) minimum TM threshold
#   TM_HIGH     (optional) maximum TM threshold
#   N_PROBES    (optional) probes per gene

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
conda activate multi_padlock_design

repo_root="$(python - <<'PY'
from multi_padlock_design import config

print(config.REPO_ROOT.resolve())
PY
)"
echo "Repository root: $repo_root"
cd "$repo_root"

FILE_PATH="$INPUT"                              # already absolute
FILE_BASE="$(basename "$FILE_PATH")"
FILE_STEM="${FILE_BASE%.*}"                     # drop .csv/.fasta

OUT_DIR="$(python - "$PARENT" "$FILE_STEM" <<'PY'
from pathlib import Path
import sys
from multi_padlock_design import config

parent, stem = sys.argv[1:3]
dest = Path(config.probe_output_root) / parent / stem
print(dest.resolve())
PY
)"
mkdir -p "${OUT_DIR}"

echo "Input file: ${FILE_PATH}"
echo "Mode: ${INPUT_TYPE}"
echo "Output dir: ${OUT_DIR}"
which python

mapfile -t PROBE_DEFAULTS < <(
python - <<'PY'
from multi_padlock_design import config

defaults = config.probe_design_defaults
print(defaults["species"])
print(defaults["arm_len"])
print(defaults["total_len"])
print(defaults["interval"])
print(defaults["tm_low"])
print(defaults["tm_high"])
print(defaults["n_probes"])
PY
)

SPECIES="${SPECIES:-${PROBE_DEFAULTS[0]}}"
ARM_LEN="${ARM_LEN:-${PROBE_DEFAULTS[1]}}"
TOTAL_LEN="${TOTAL_LEN:-${PROBE_DEFAULTS[2]}}"
INTERVAL="${INTERVAL:-${PROBE_DEFAULTS[3]}}"
TM_LOW="${TM_LOW:-${PROBE_DEFAULTS[4]}}"
TM_HIGH="${TM_HIGH:-${PROBE_DEFAULTS[5]}}"
N_PROBES="${N_PROBES:-${PROBE_DEFAULTS[6]}}"

# Feed probedesign.py prompts
if [[ "${INPUT_TYPE,,}" == "csv" ]]; then
  python multi_padlock_design/probedesign.py <<- ENDOF
${OUT_DIR}
${SPECIES}
${INPUT_TYPE,,}
${FILE_PATH}
${ARM_LEN}
${TOTAL_LEN}
${INTERVAL}
${TM_LOW}
${TM_HIGH}
${N_PROBES}
ENDOF
elif [[ "${INPUT_TYPE,,}" == "fasta" ]]; then
  python multi_padlock_design/probedesign.py <<- ENDOF
${OUT_DIR}
${SPECIES}
${INPUT_TYPE,,}

${FILE_PATH}
${ARM_LEN}
${TOTAL_LEN}
${INTERVAL}
${TM_LOW}
${TM_HIGH}
${N_PROBES}
ENDOF
else
  echo "Unknown INPUT_TYPE: ${INPUT_TYPE}" >&2
  exit 1
fi
