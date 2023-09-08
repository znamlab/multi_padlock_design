# multi_padlock_design
Padlock design tool modified for running on slurm

First create a .csv containing your genes of interest. <br />
Run fasta_generator.py inputting the .csv file, to generate a single fasta for genes of interest. <br />
Split this single fasta file into a set of individual fasta files that are grouped by gene using split.sh <br />
Then edit probe_design and probe_design_parallel.sh and run parallel in order to create a set of jobs for each gene. Each will take a few hours to run. <br />
Collate the outputs for each gene from the scratch dir using collate_probes <br />
Substitute the N10 barcode region for barcodes selected at random from a barcode.txt file using select_barcodes <br />
Delete scratch files with autoremovescratchpadlocks <br />
