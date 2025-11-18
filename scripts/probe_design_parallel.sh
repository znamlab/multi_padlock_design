#!/bin/bash
#SBATCH --job-name=probe_design_parallel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --partition=ncpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.becalick@crick.ac.uk

set -euo pipefail

# Usage: probe_design_parallel.sh <input_file.(fa|fasta|csv)>
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input_file.(fa|fasta|csv)>"
  exit 1
fi

if command -v ml &>/dev/null; then
  ml purge
  ml Anaconda3
else
  echo "Module system is not available on node: $CURRENT_NODE"
fi

eval "$(/camp/apps/eb/software/Anaconda3/2024.10-1/bin/conda shell.bash hook)"
conda activate multi_padlock_design

# Determine repository root (prefer git, fallback to this script's directory)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PARENT_DIR="$(cd "$SCRIPT_DIR/.." && pwd -P)"
if repo_root=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel 2>/dev/null); then
  :
elif repo_root=$(git -C "$PARENT_DIR" rev-parse --show-toplevel 2>/dev/null); then
  :
else
  repo_root="$PARENT_DIR"
fi

# Ensure we run from repo root so Python can import local config
cd "$repo_root"

# Canonicalize input file to absolute path
if command -v realpath >/dev/null 2>&1; then
  input_file="$(realpath -m "$1")"
else
  input_file="$(cd "$(dirname "$1")" && pwd -P)/$(basename "$1")"
fi

if [[ ! -f "$input_file" ]]; then
  echo "Error: input_file '$input_file' does not exist"
  exit 1
fi

ext_lower=$(echo "${input_file##*.}" | tr '[:upper:]' '[:lower:]')
case "$ext_lower" in
  fa|fasta)
    input_type="fasta"
    ;;
  csv)
    input_type="csv"
    ;;
  *)
    echo "Error: input file extension must be .fa, .fasta, or .csv"
    exit 1
    ;;
esac

echo "Preparing to split input file for processing"
# Split input into per-item files in scratch
if [[ "$input_type" == "fasta" ]]; then
  python "${repo_root}/multi_padlock_design/io/fasta_splitter.py" "$input_file"
else
  python "${repo_root}/multi_padlock_design/io/csv_splitter.py" "$input_file"
fi

# Compute output directory based on Python config.split_input and input stem
out_dir="$(python - "$input_file" <<'PY'
from pathlib import Path
import sys
from multi_padlock_design import config
inp = Path(sys.argv[1])
print((config.split_input / inp.stem).resolve())
PY
)"

parent_base="$(basename "$out_dir")"
logdir="${repo_root}/logs/slurm_logs/${parent_base}"
mkdir -p "$logdir"

# Iterate over produced files in scratch directory and submit jobs
shopt -s nullglob nocaseglob
if [[ "$input_type" == "fasta" ]]; then
  files=( "${out_dir}"/*.fasta )
else
  files=( "${out_dir}"/*.csv )
fi
shopt -u nocaseglob

if (( ${#files[@]} == 0 )); then
  echo "No split files found in ${out_dir}"
  exit 0
fi

for src in "${files[@]}"; do
  base_noext="$(basename "${src%.*}")"
  echo "${parent_base}"
  echo "${base_noext}"
  echo "Starting job ${src}"
  sbatch \
    --export=INPUT="$src",PARENT="$parent_base",INPUT_TYPE="$input_type" \
  --output="${logdir}/${parent_base}_${base_noext}.out" \
  "${repo_root}/scripts/probe_design.sh"
done
