# multi_padlock_design
Padlock design tool modified for running on slurm

First create a .csv containing your genes of interest. <br />
Run csv_splitter.py inputting the .csv file, to generate a single csv for each of the genes of interest. <br />
Then run probe_design_parallel.sh, inputting the directory containing the gene csvs as the only argument. This will run parallel jobs for each gene. Each will take a few hours to run. <br />
Collate the outputs for each gene from the scratch dir using collate_probes <br />
Substitute the N10 barcode region for barcodes selected at random from a barcode.txt file using select_barcodes <br />
Delete scratch files with autoremovescratchpadlocks <br />
