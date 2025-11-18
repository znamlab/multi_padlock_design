# Multi Padlock Design Pipeline

A padlock design workflow for the generation of in situ sequencing padlock probes.
Based originally on code developed in the Mats Nilsson Lab, Stockholm University (Xiaoyan, 2018), expanded to include support for SLURM job parallelisation, the use of Ensembl BLAST databases, non CDS regions, melting temperature based specificity criteria, barcode and hybridisation region design, secondary structure analysis by the lab of Petr Znamenskiy (Becalick, 2025).

## Installation and Setup

First clone the repository and install it with with conda (this will automatically install all required packages):

```bash
git clone git@github.com:znamlab/multi_padlock_design.git
cd multi_padlock_design
conda env create -f environment.yml
```
Then, activate the new environment:

```bash
conda activate multi-padlock-design
```

Next, install the following:
- Install BLAST+ and add it to your `PATH`.
- Install ClustalW2 and add it to your `PATH`.
- Install Clustal-Omega and add it to your `PATH`.


Finally, download the reference databases for BLAST searching. You can do this automatically by running the following command once you've activated your conda environment:
```bash
mpd-download-databases
```
Use `mpd-download-databases --help` to explore options for download.
If a checksum mismatch occurs, this is likely due to an interruption during download from the FTP server.
If you want to select a different version of your refseq/ensembl database or use a different species, you can do so by downloading the data manually from the following locations:
- Download RefSeq mouse/human mRNA sequences from `https://ftp.ncbi.nlm.nih.gov/refseq/`
- You need all files from mRNA_Prot ending with .rna.fna.gz
- Download Ensembl mouse/human CDS and cDNA sequences from `https://ftp.ensembl.org/pub/release-115/fasta/`
- You need the files ending in .all.fa.gz from the cds and cdna subfolders
- If using Ensembl, also download the gencode annotations from https://ftp.ebi.ac.uk/pub/databases/gencode/
- You need the gencode.v*.basic.annotation.gtf.gz file from the respective species' latest release
- Unzip all the .gz files
- Place these files into subfolders called `refseq` and `ensembl`

## Configuration
- Update `multi_padlock_design/config.py` with the paths to the downloaded mRNA sequence files.
- Set the maximum number of parallel threads for multiple sequence alignment (MSA) and BLAST searches in `config.py`.
- You should also select your reference transcriptome of choice. Ensembl is required if you want to design probes against non CDS regions and blast all probes against these regions for off-target searching.

## Usage
If using SLURM:

If running locally:

## Supported Input Files
1. CSV file containing gene acronyms with the corresponding linker and barcode sequences, one gene per row, no header.
2. FASTA file containing target sequences, each with their own FASTA header.

## Adjustable Parameters
- **Species**: Reference database used for specificity checks and for resolving gene acronyms when FASTA inputs are absent.
- **Padlock arm length**: Length of each target arm (final target sequence length is twice this value).
- **Tm range**: Lower and upper melting temperature bounds, assuming 0.1 µM probe, 0.075 M monovalent salt, 0.01 M bivalent salt, and 20% formamide.

## Specificity Criteria
A target sequence is considered sufficiently specific only if no sequence in the reference database exceeds:
- A melting temperature of 37C for both the left **and** right padlock arm regions (assumed to be necessary for a ligation event to occur)
- No gaps or mismatches within 3bp of the ligation site.
