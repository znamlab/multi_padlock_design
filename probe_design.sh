#!/bin/bash
#SBATCH --job-name=probe_design
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --partition=ncpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.becalick@crick.ac.uk
#SBATCH --output=/nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/${PARENT}/${PARENT}_%j.out

# Get the current node name
CURRENT_NODE=$(hostname)

# Function to resubmit job excluding the current node
resubmit_job() {
    echo "Resubmitting job excluding node: $CURRENT_NODE"
    sbatch --exclude=$CURRENT_NODE "$0"
    exit 0
}

# Check if the 'ml' command is available
if command -v ml &> /dev/null
then
    ml purge
    ml Anaconda3/2024.10-1
    ml Clustal-Omega
    ml BLAST+
    ml ClustalW2
else
    echo "Module system is not available on node: $CURRENT_NODE"
    resubmit_job
fi

# Check if symbolic link target directory exists before creating the link
if [ -d "/nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/${PARENT}/" ]; then
  ln -f /nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/${PARENT}/${PARENT}_${SLURM_JOB_ID}.out /nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/${PARENT}/${PARENT}_${INPUT}.out
fi

eval "$(/camp/apps/eb/software/Anaconda3/2024.10-1/bin/conda shell.bash hook)"
conda activate iss-preprocess
cd /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design

# Printing the file path
echo /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/${PARENT}/${INPUT}
which python
# Running the Python script with the specified parameters
python probedesign.py <<- ENDOF
mouse
/nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/${PARENT}/${INPUT}
/nemo/lab/znamenskiyp/scratch/${INPUT}
20
8
60
78
20
ENDOF

# Remove symbolic link after job completion
if [ -L "/nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/${PARENT}/${PARENT}_%j.out" ]; then
  rm /nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/${PARENT}/${PARENT}_%j.out
fi
