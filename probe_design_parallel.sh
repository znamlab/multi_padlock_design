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

# Usage: probe_design_parallel.sh <input_dir> <fasta|csv>
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input_dir> <fasta|csv>"
  exit 1
fi

# Canonicalize input_dir to an absolute path (portable fallback if realpath is missing)
if command -v realpath >/dev/null 2>&1; then
  input_dir="$(realpath -m "$1")"
else
  # -P resolves symlinks; cd handles relative paths
  input_dir="$(cd "$1" && pwd -P)"
fi
input_type="$(echo "$2" | tr '[:upper:]' '[:lower:]')"

if [[ ! -d "$input_dir" ]]; then
  echo "Error: input_dir '$input_dir' does not exist"
  exit 1
fi
if [[ "$input_type" != "fasta" && "$input_type" != "csv" ]]; then
  echo "Error: input_type must be 'fasta' or 'csv'"
  exit 1
fi

parent_base="$(basename "$input_dir")"
logdir="/nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/${parent_base}"
mkdir -p "$logdir"

# Find files directly in input_dir (no copying)
shopt -s nullglob nocaseglob
files=( "${input_dir}"/*.${input_type} )
shopt -u nocaseglob

if (( ${#files[@]} == 0 )); then
  echo "No .${input_type} files found in ${input_dir}"
  exit 0
fi

for src in "${files[@]}"; do
  base_noext="$(basename "${src%.*}")"
  echo "Starting job ${src}"
  sbatch \
    --export=INPUT="$src",PARENT="$parent_base",INPUT_TYPE="$input_type" \
    --output="${logdir}/${parent_base}_${base_noext}.out" \
    /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/probe_design.sh
done
