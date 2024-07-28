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

#tmux
sleep 3

ml Anaconda3/2022.05
ml Clustal-Omega
ml BLAST+
ml ClustalW2

# Creating symbolic link for the output log
ln -f /nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/${PARENT}/${PARENT}_%j.out /nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/${PARENT}/${PARENT}_${INPUT}.out

source /camp/apps/eb/software/Anaconda/conda.env.sh
conda activate base
cd /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design

# Printing the file path
echo /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/${PARENT}/${INPUT}

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

# Removing the symbolic link after the job is done
rm /nemo/lab/znamenskiyp/home/users/becalia/logs/slurm_logs/${PARENT}/${PARENT}_%j.out

