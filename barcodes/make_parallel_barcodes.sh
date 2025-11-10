#!/bin/bash
#SBATCH --job-name=make_barcodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00
#SBATCH --mem=8G
#SBATCH --partition=ncpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=colasa@crick.ac.uk

ml purge
ml Anaconda3/2024.10-1
conda activate spapros

# prevent oversubscription (each worker is a process)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd /nemo/lab/znamenskiyp/home/users/colasa/code/multi_padlock_design

python make_nested_barcodes_parallel.py \
  --input /nemo/lab/znamenskiyp/home/users/colasa/code/multi_padlock_design/constrained_panel/1800_d4_10bp_barcodes_improved.txt \
  --outdir constrained_panel \
  --processes 16 \
  --base-seed 50 \
  --iters 8000 \
  --subset-restarts 8 \
  --alpha 3.0 \
  --d-global 4 \
  --d-prefix 3
