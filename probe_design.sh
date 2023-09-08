#!/bin/bash
#SBATCH --job-name=probe_design
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --partition=cpu
#SBATCH --mail-type=ALL
#SBATCH -o /nemo/lab/znamenskiyp/home/users/becalia/logs/probe_design.out
#SBATCH -e /nemo/lab/znamenskiyp/home/users/becalia/logs/probe_design.err
#SBATCH --mail-user=alexander.becalick@crick.ac.uk


#tmux
sleep 3

ml Anaconda3/2022.05
ml Clustal-Omega
ml BLAST+
ml ClustalW2

source /camp/apps/eb/software/Anaconda/conda.env.sh
conda activate base
cd /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design

python probedesign.py << ENDOF
mouse

/nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/taar_genes/$INPUT
/nemo/lab/znamenskiyp/scratch/$INPUT
20
4
60
78
20
ENDOF