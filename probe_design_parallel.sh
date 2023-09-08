#!/bin/bash
#SBATCH --job-name=probe_design_parallel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --partition=cpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.becalick@crick.ac.uk


#tmux
sleep 3
cd /nemo/lab/znamenskiyp/home/users/becalia/code/multi_padlock_design/taar_genes/

for f in *.fasta
do
 	echo Starting job "$f"
	sbatch --export=INPUT="$f" /nemo/lab/znamenskiyp/home/users/becalia/code/SBATCH/probe_design.sh --output=/nemo/lab/znamenskiyp/home/users/becalia/logs/taar_genes_"$f".out
done
