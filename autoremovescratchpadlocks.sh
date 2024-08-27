#!/bin/bash
#SBATCH --job-name=probe_design
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=8G
#SBATCH --partition=ncpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.becalick@crick.ac.uk


# Navigate to the root directory
cd /nemo/lab/znamenskiyp/scratch/

# Find '*.csv' directories
for dir in *.csv/; do
  if [ -d "$dir" ]; then
    # Find and delete 'TempFolder*' directories within  directories
    find "$dir" -type d -name 'TempFolder*' -exec rm -r '{}' \;
  fi
done
