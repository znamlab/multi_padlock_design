#!/bin/bash
#SBATCH --job-name=probe_design_parallel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --partition=ncpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.becalick@crick.ac.uk

# Check if the input parameter is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <input>"
  exit 1
fi

csv_folder=$1

#tmux
sleep 3
mkdir -p /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/${csv_folder}/
cd /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/${csv_folder}/

for file in *.csv
do
  echo "Starting job $file"
  sbatch --export=INPUT="$file",PARENT="$csv_folder" --output=/nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/${csv_folder}/${csv_folder}_"${file%.csv}".out /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/probe_design.sh
done
