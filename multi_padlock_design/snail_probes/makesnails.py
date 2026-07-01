from pathlib import Path
import pandas as pd
import numpy as np
import random
import itertools
import lib.readblast as readblast  
import re
import matplotlib.pyplot as plt



def revcomp(seq):
    """Return the reverse complement of a DNA sequence string (A/T/C/G only)."""
    complement = str.maketrans('ATCG', 'TAGC')
    return seq.translate(complement)[::-1]

def read_candidates(outdir, panel_name, probefile = "3.AllSpecificTargets_*.csv"):
    """
    Read and concatenate all "3.AllSpecificTargets_*.csv" files from subdirectories of the given panel directory.
    
    Args:
        outdir (Path): Base output directory containing the panel subdirectory.
        panel_name (str): Name of the panel subdirectory.

    Returns:
        pd.DataFrame: Concatenated dataframe with columns ['target', 'Tm', 'startpos', 'gene'].
    """
    all_dfs = []
    panel_dir = outdir / panel_name

    # Loop over subdirectories
    for subdir in panel_dir.iterdir():
        if subdir.is_dir():
            print(f"Processing {subdir.name}")

            # Collect "3.AllSpecificTargets_*.csv"
            target_files = list(subdir.glob(probefile))
            if len(target_files) != 1:
                print(f"Warning: Expected exactly one target file in {subdir}, found {len(target_files)}. Skipping.")
                continue
            tf = target_files[0]

            # Load and tag with gene name (= subdir name)
            with open(tf) as f:
                clean_lines = [line for line in f if not line.startswith(">")]
                clean_lines = [smart_split(line.strip()) for line in clean_lines]

            header, *data = clean_lines
            if len(data)==0:
                print(f"Warning: No data found in {tf}, skipping.")
                continue
            if len(data[0])<2:
                print(f"Warning: No data found in {tf}, skipping.")
                continue
            df = pd.DataFrame(data, columns=header)
            df["gene"] = subdir.name
            if len(df)<1:
                print(f"Warning: No targets found in {subdir.name}, skipping.")
                continue
            all_dfs.append(df)

    # Concatenate all into one big dataframe (OUTSIDE the loop)
    if all_dfs:
        big_df = pd.concat(all_dfs, ignore_index=True)
    else:
        big_df = pd.DataFrame()  # Return empty DataFrame if no data
    
    big_df = big_df[big_df['acronym'].notna() & (big_df['acronym'] != '')]

    return big_df

def smart_split(line):
    """Split a csv line on commas, ignoring commas nested inside [ ... ] fields."""
    return re.split(r',(?![^[]*\])', line)

rows = []


def generate_random_dna_sequences(min_barcode, length, num_sequences, for_overlap, genes=None):
    """
    Randomly generate a set of mutually distinguishable DNA sequences, used as the
    shared overlap region between a padlock and its matching primer.

    Sequences are rejected and re-drawn until they satisfy all of:
      - Tm (nearest-neighbour) between 25 and 35 C.
      - No homopolymer run of 3+ identical bases.
      - If for_overlap is True, bases [3:6] avoid ligation-blocking motifs
        (their reverse complement would present a TXA at the ligation junction).
      - At least 2 mismatches against every other sequence's first min_barcode
        bases, and at least 4 mismatches against every other full sequence
        (keeps sequences distinguishable even with sequencing errors).

    Args:
        min_barcode (int): Number of leading bases used for the stricter
            (>=2 mismatch) distinguishability check.
        length (int): Length of each generated sequence.
        num_sequences (int): How many distinct sequences to generate.
        for_overlap (bool): Whether to apply the ligation-junction filter.
        genes (list[str], optional): Gene/padlock labels to attach to each row
            (padded with "no_gene" or truncated to match the number of rows).

    Returns:
        pd.DataFrame: One row per sequence, columns ['Sequence',
        'Reverse Complement', 'Tms', 'Overlap Primer Sequence', 'ID'].
    """

    bases = ['A', 'T', 'C', 'G']
    sequences = []
    reverse_complements = []
    Tms = []
    
    while len(sequences) < num_sequences:
        sequence = ''.join(random.choice(bases) for _ in range(length))
        
        # Check condition for overlap oligos
        if for_overlap and sequence[3:6] in ['TAA', 'TCA', 'TGA', 'TTA', 'ATT', 'AGT', 'ACT', "AAT"]:#ligation wont happen when reverse complement is TXA In 6:9 position for overlap
            continue
        
        #No homopolymers
        if ('AAA' in sequence) or ('TTT' in sequence) or ('CCC' in sequence) or ('GGG' in sequence):
            continue
            
        # check Tm?
        Tm = readblast.calc_tm_NN(sequence)
        if Tm < 25 or Tm > 35:
            continue
        
        if all(count_diff(seq[:min_barcode], sequence[:min_barcode]) >= 2 for seq in sequences) and \
           all(count_diff(seq, sequence) >= 4 for seq in sequences):
            sequences.append(sequence)
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            reverse_seq = sequence[::-1]
            reverse_complement_seq = ''.join(complement[base] for base in reverse_seq)
            reverse_complements.append(reverse_complement_seq)
            Tms.append(Tm)

    # Build dataframe
    df = pd.DataFrame({
        'Sequence': sequences,
        'Reverse Complement': reverse_complements, 
        'Tms': Tms
    }).drop_duplicates()
    
    # Handle gene names
    if genes is None:
        genes = ["no_gene"] * len(df)
    else:
        # pad or truncate to match df length
        if len(genes) < len(df):
            genes = genes + ["no_gene"] * (len(df) - len(genes))
        elif len(genes) > len(df):
            genes = genes[:len(df)]
    
    df['Overlap Primer Sequence'] = genes
    df['ID'] = range(len(df))  # unique IDs
    
    return df

def count_diff(seq1, seq2):
    """Count the number of mismatched positions between two equal-length sequences."""
    return sum(base1 != base2 for base1, base2 in zip(seq1, seq2))

def get_barcodes(csv_path, gene_number=None, indices=None):
    """
    Read a barcodes csv (columns: barcode, is_7mer, gene, panel).

    If indices is given, return that subset of rows.
    Otherwise, return the first gene_number unused (gene is empty) rows.
    """
    df = pd.read_csv(csv_path)

    if indices is not None:
        return df.loc[indices]

    unused = df[df['gene'].isna() | (df['gene'] == '')]
    return unused.iloc[:gene_number]


class snailmaker():
    """
    Builds a full SNAIL (padlock + primer) probe panel for a set of genes.

    Given candidate target sequences (from read_candidates), this class picks
    the best targets per gene, generates unique overlap sequences linking each
    padlock to its primer, assigns barcodes/UMIs for sequencing-based readout
    (or a fixed hybridization probe for hyb_genes), and assembles the final
    full-length oligo sequence for every padlock and primer. See
    docs/probe_construction.md for a full description of the pipeline and the
    resulting probe layout.

    Call main() to run the whole pipeline end-to-end, or call the individual
    methods (assign_targets, generate_extra_redundancy, assemble_probes, ...)
    directly for more control.
    """

    def __init__(self,
                 outdir,
                 panel_name,
                 hyb_genes,
                 hyb_probes_path,
                 barcode_csv_path,
                 probefile = "3.AllSpecificTargets_*.csv",
                 plot_report = True,
                 min_probes = 4,
                 sequencing_primer = 'GTCGGACTGTAGAACTCT'
                 ):
        """
        Args:
            outdir (Path): Base output directory containing the panel subdirectory.
            panel_name (str): Name of the panel subdirectory (holds one
                subdirectory of BLAST output per gene).
            hyb_genes (list[str]): Genes to be read out by a fixed hybridization
                probe (from hyb_probes_path) rather than a sequencing barcode.
            hyb_probes_path (str | Path): CSV of available hybridization probes
                (columns include 'SNAIL TARG' and 'SEQ').
            barcode_csv_path (str | Path): CSV of available sequencing barcodes
                (columns: barcode, is_7mer, gene, panel) to draw from for the
                non-hyb genes.
            probefile (str): Glob pattern for the per-gene BLAST target file.
            plot_report (bool): If True, generated figures (e.g. probes_per_gene)
                are displayed when show_figures() is called.
            min_probes (int): Number of padlock/primer pairs per gene (max 4,
                since each pair is tagged with a single-base UMI).
            sequencing_primer (str): Constant sequencing primer sequence baked
                into every non-hyb padlock.
        """
        if min_probes > 4:
            raise ValueError("min_probes cannot exceed 4: probes are UMI-tagged with a single base (A/T/C/G)")

        self.hyb_genes = hyb_genes
        self.plot_report = plot_report
        self.min_probes = min_probes
        self.sequencing_primer = sequencing_primer
        self.big_df = self.read_candidates(outdir, panel_name, probefile)
        self.hyb_probes = pd.read_csv(hyb_probes_path)
        self.barcode_csv_path = barcode_csv_path
        self.figures = []
        self.dropped_genes = []
        self.outdir = Path(outdir) / panel_name
        self.panel_name = panel_name

    def read_candidates(self, outdir, panel_name, probefile = "3.AllSpecificTargets_*.csv"):
        """Instance wrapper around the module-level read_candidates (see its docstring)."""
        return read_candidates(outdir, panel_name, probefile)

    def probes_per_gene(self):
        """
        Plot the number of candidate targets found per gene (self.big_df),
        with a horizontal line at min_probes marking the minimum needed to
        keep a gene in the panel. Stores the figure in self.figures and
        returns it.
        """
        per_gene = self.big_df.groupby('acronym')['target'].count()

        # Sort values for clarity
        per_gene_sorted = per_gene.sort_values(ascending=False)

        # Make figure tall enough for ~100 bars
        fig, ax = plt.subplots(figsize=(len(per_gene_sorted) * 0.3, 10))

        per_gene_sorted.plot(kind='bar', ax=ax)

        ax.set_xlabel("Gene Acronym")
        ax.axhline(y = self.min_probes, color = 'red')
        ax.set_ylabel("Target Count")
        ax.set_title("Target Counts per Gene Acronym")
        fig.tight_layout()

        self.figures.append(fig)
        return fig

    def show_figures(self):
        """Display every figure collected in self.figures, but only if plot_report is True."""
        if not self.plot_report:
            return

        for fig in self.figures:
            plt.figure(fig.number)
            plt.show()

    def choose_targets(self, gene, n_probes=None):
        """
        Pick the best n_probes targets within a single gene, preferring
        targets with fewer/shorter homopolymer runs.
        """
        if n_probes is None:
            n_probes = self.min_probes

        candidates = self.big_df[self.big_df['acronym'] == gene].copy()
        candidates['homopolymer_score'] = np.zeros(len(candidates))

        for idx, row in candidates.iterrows():
            target = row['target']
            if ('AAAAA' in target) or ('TTTTT' in target) or ('CCCCC' in target) or ('GGGGG' in target):
                candidates.at[idx, 'homopolymer_score'] = 0 #never length 5 homopolymers if you can avoid it
            elif ('AAAA' in target) or ('TTTT' in target) or ('CCCC' in target) or ('GGGG' in target):
                candidates.at[idx, 'homopolymer_score'] = 1 #avoid length 4 homopolymers if you can
            else:
                candidates.at[idx, 'homopolymer_score'] = 2 #ideal, no homopolymers of length 4 or greater

        candidates.sort_values(by=['homopolymer_score'], ascending=False, inplace=True)
        chosen = candidates.iloc[0:n_probes]
        return chosen

    def build_overlaps(self):
        """
        Expand self.overlaps (one row per padlock) into one row per padlock
        AND one matching row per primer, so each padlock/primer pair shares
        the same overlap sequence. Primer rows are the padlock's row with an
        'Overlap Primer Sequence' ID offset by +10 (e.g. GENE_1 -> GENE_11),
        matching the padlock/primer numbering used elsewhere (see main()).
        """
        # Create a list to store all rows (original + duplicates)
        new_rows = []

        # Iterate through the dataframe in groups of min_probes
        for i in range(0, len(self.overlaps), self.min_probes):
            # Get the current group of min_probes rows
            group = self.overlaps.iloc[i:i+self.min_probes]
            
            # Add the original 4 rows
            new_rows.append(group)
            
            # Create duplicates with modified IDs
            duplicates = group.copy()
            duplicates['Overlap Primer Sequence'] = duplicates['Overlap Primer Sequence'].str.replace(
                r'_(\d+)$', 
                lambda m: f'_{int(m.group(1)) + 10}',
                regex=True
            )
            
            # Add the duplicate 4 rows
            new_rows.append(duplicates)

        # Concatenate all rows and reset index
        self.overlaps = pd.concat(new_rows, ignore_index=True)

    def assign_targets(self):
        """
        For every gene in self.genes, pick min_probes targets (via
        choose_targets) and write their target/Tm/startpos/endpos/
        transcript_region into the matching padlock and primer rows of
        self.probes_df. Genes with zero candidates, or fewer candidates than
        min_probes, are dropped entirely (removed from self.genes and
        self.probes_df, recorded in self.dropped_genes) rather than being
        assigned a partial/undersized set of probes.
        """
        for gene in list(self.genes):
            candidates = self.big_df[self.big_df['acronym'] == gene]
            if candidates.empty:
                print(f"Warning: no targets for {gene}, gene DROPPED")
                self.probes_df = self.probes_df[self.probes_df['gene'] != gene].copy()
                self.genes = [g for g in self.genes if g != gene]
                self.dropped_genes.append(gene)
                continue
            if candidates.shape[0] < self.min_probes:
                print(f"Warning: only {candidates.shape[0]} targets for {gene}, fewer than {self.min_probes}, gene DROPPED")
                self.probes_df = self.probes_df[self.probes_df['gene'] != gene].copy()
                self.genes = [g for g in self.genes if g != gene]
                self.dropped_genes.append(gene)
                continue

            n = min(self.min_probes, len(candidates))
            chosen = self.choose_targets(gene, n_probes=n)

            for pair in range(1, self.min_probes + 1):
                padlock = f"{gene}_{pair}"
                primer = f"{gene}_{pair + 10}"
                target = chosen.iloc[pair - 1]

                for probe_id, probe_type in [(padlock, 'padlock'), (primer, 'primer')]:
                    self.probes_df.loc[probe_id, 'target'] = target['target']
                    self.probes_df.loc[probe_id, 'Tm'] = target['Tm']
                    self.probes_df.loc[probe_id, 'startpos'] = target['startpos']
                    self.probes_df.loc[probe_id, 'endpos'] = target['endpos']
                    self.probes_df.loc[probe_id, 'transcript_region'] = target['transcript_region']
                    self.probes_df.loc[probe_id, 'type'] = probe_type

    def generate_extra_redundancy(self):
        """
        Assign each gene a random 5-mer barcode (no length-3 homopolymers)
        for extra redundancy.
        """
        all_combinations = [''.join(combo) for combo in itertools.product('ATCG', repeat=5)]

        homopolymers = ['AAA', 'TTT', 'CCC', 'GGG']
        filtered_combinations = [seq for seq in all_combinations if not any(hp in seq for hp in homopolymers)]

        random_selection = random.sample(filtered_combinations, len(self.genes))

        self.extra_redundancy = pd.DataFrame({
            'genes': self.genes,
            'barcode': random_selection
        }).set_index('genes')

        return self.extra_redundancy

    def assemble_probes(self):
        """
        Build the final full-length oligo sequence for every padlock/primer pair.

        Padlock layout (5' -> 3'):
          hyb readout:        overlap_5' + target_seq + hyb_sequence + extra_redundancy[0:2] + overlap_3'
          sequencing readout: overlap_5' + target_seq + sequencing_primer + barcode + UMI + extra_redundancy + overlap_3'

        Primer layout (5' -> 3'), same for both readout types:
          target_seq + overlap_sequence

        Where:
          - overlap_5'/overlap_3' are the two halves of the shared padlock/primer
            overlap (from build_overlaps), reverse-complemented so the padlock's
            ends are complementary to the target flanking the primer's binding site.
          - target_seq is the reverse complement of the first 22 bases of the
            target (padlock) or of bases [23:] (primer) — the 1-base gap between
            them is intentional (see the "leave a distance of 1" comment).
          - hyb_sequence is the padlock-specific hybridization readout sequence
            (only for genes in hyb_genes).
          - barcode/UMI/extra_redundancy identify the gene and probe pair for
            genes read out by sequencing.

        Requires self.genes, self.probes_df (with 'hyb', 'target',
        'overlap_sequence' columns), self.extra_redundancy, self.gene_dict, and
        self.hyb_probes to already be populated (see main()).
        """
        self.probes_df['sequencing_primer'] = self.sequencing_primer

        new_probes = []
        UMI = ['A', 'T', 'C', 'G']
        for gene in self.genes:
            subset_probes = self.probes_df.loc[(self.probes_df['gene'] == gene).values].copy()
            for pair in range(1, self.min_probes + 1):
                padlock = f"{gene}_{pair}"
                primer = f'{gene}_{pair+10}'
                print(f'Assembling {padlock}')

                #padlock bits
                subset_probes.loc[padlock, 'extra_redundancy'] = self.extra_redundancy.loc[gene].values[0]
                subset_probes.loc[padlock, 'target_seq'] = revcomp(subset_probes.loc[padlock, 'target'][0:22])

                if subset_probes.loc[padlock, 'hyb']: #readout by hyb
                    subset_probes.loc[padlock, 'hyb_sequence'] = revcomp(self.hyb_probes[self.hyb_probes['SNAIL TARG']==gene]['SEQ'].values[0].split('/')[2])
                else: #readout by seq
                    subset_probes.loc[padlock, 'barcode'] = self.gene_dict[self.gene_dict['genes']==gene]['barcodes'].values[0]
                    subset_probes.loc[padlock, 'is_7_mer'] = self.gene_dict[self.gene_dict['genes']==gene]['is_7_mer'].values[0]
                    subset_probes.loc[padlock, 'UMI'] = UMI[pair-1]

                subset_probes.loc[padlock, 'overlap_5_seq'] =  revcomp(subset_probes.loc[padlock, 'overlap_sequence'][0:6])
                subset_probes.loc[padlock, 'overlap_3_seq'] =  revcomp(subset_probes.loc[padlock, 'overlap_sequence'][6:])

                #primer bits
                subset_probes.loc[primer, 'target_seq'] = revcomp(subset_probes.loc[primer, 'target'][23:]) #leave a distance of 1 in the target

                #now, assemble
                def access(field, kind = 'padlock', reference = {'padlock':padlock, 'primer':primer}, subset_probes = subset_probes):
                    id = reference[kind]
                    return str(subset_probes.loc[id, field])

                if subset_probes.loc[padlock, 'hyb']: #read out by hybridization
                    pad_sequence = access('overlap_5_seq') + access('target_seq') + access('hyb_sequence') + access('extra_redundancy')[0:2] + access('overlap_3_seq') #ensure RCA optimal length
                    subset_probes.loc[padlock, 'full_sequence'] = pad_sequence
                    subset_probes.loc[padlock, 'length'] = len(pad_sequence)

                    prim_sequence = access('target_seq', 'primer') + access('overlap_sequence', 'primer')
                    subset_probes.loc[primer, 'full_sequence'] = prim_sequence
                    subset_probes.loc[primer, 'length'] = len(prim_sequence)

                    continue
                else: #read out by sequencing

                    #padlocks
                    pad_sequence = access('overlap_5_seq') + access('target_seq') + access('sequencing_primer') + access('barcode') + access('UMI') + access('extra_redundancy') + access('overlap_3_seq')
                    subset_probes.loc[padlock, 'full_sequence'] = pad_sequence
                    subset_probes.loc[padlock, 'length'] = len(pad_sequence)

                    #primers
                    prim_sequence = access('target_seq', 'primer') + access('overlap_sequence', 'primer')
                    subset_probes.loc[primer, 'full_sequence'] = prim_sequence
                    subset_probes.loc[primer, 'length'] = len(prim_sequence)

            new_probes.append(subset_probes)

        self.probes_df = pd.concat(new_probes)
        return self.probes_df

    def save_probes(self, outdir):
        """Write self.probes_df (the final assembled probe set) to outdir/probes_df.csv."""
        panel_dir = Path(outdir)
        panel_dir.mkdir(parents=True, exist_ok=True)
        out_path = panel_dir / "probes_df.csv"
        self.probes_df.to_csv(out_path)
        return out_path

    def save_assigned_barcodes(self, outdir):
        """
        Write self.gene_dict (gene -> barcode/is_7_mer assignment for this
        panel) to outdir/assigned_barcodes.csv. This does NOT write back into
        the source barcode_csv_path — barcode usage is only tracked per-panel,
        here.
        """
        panel_dir = Path(outdir)
        panel_dir.mkdir(parents=True, exist_ok=True)
        out_path = panel_dir / "assigned_barcodes.csv"
        self.gene_dict.to_csv(out_path)
        return out_path

    def write_summary(self, outdir):
        """Write a human-readable summary.txt (gene/barcode/dropped-gene counts) to outdir."""
        panel_dir = Path(outdir)
        panel_dir.mkdir(parents=True, exist_ok=True)
        out_path = panel_dir / "summary.txt"

        n_barcodes = len(self.gene_dict)
        n_7mers = int(self.gene_dict['is_7_mer'].sum())

        lines = [
            f"Number of genes: {len(self.genes)}",
            f"Genes: {', '.join(self.genes)}",
            f"Number of hyb probes: {len(self.hyb_probes)}",
            f"Genes with hyb probes: {', '.join(self.hyb_genes)}",
            f"Number of barcodes: {n_barcodes}",
            f"Number of 7-mers within the barcodes: {n_7mers}",
            f"Number of dropped genes: {len(self.dropped_genes)}",
            f"DROPPED GENES: {', '.join(self.dropped_genes)}",
        ]

        with open(out_path, 'w') as f:
            f.write('\n'.join(lines) + '\n')

        return out_path

    def main(self):
        """
        Run the full probe-design pipeline end-to-end:
          1. Plot targets-per-gene (if plot_report).
          2. Generate unique padlock/primer overlap sequences (build_overlaps).
          3. Build the empty probes_df skeleton (one row per padlock/primer).
          4. Assign the best target per padlock/primer (assign_targets),
             dropping genes with insufficient candidates.
          5. Assign sequencing barcodes to non-hyb genes (gene_dict).
          6. Generate per-gene extra-redundancy barcodes.
          7. Assemble the final oligo sequences (assemble_probes).
          8. Save probes_df.csv, assigned_barcodes.csv, and summary.txt to
             self.outdir, then display any collected figures.

        Returns:
            pd.DataFrame: The final assembled self.probes_df.
        """
        if self.plot_report:
            self.probes_per_gene()

        self.genes = self.big_df['acronym'].unique().tolist()

        padlock_ids = list(range(1, self.min_probes + 1))
        primer_ids = [pid + 10 for pid in padlock_ids]
        probe_ids = padlock_ids + primer_ids

        probes = np.array([f'{gene}_{pid}' for gene in self.genes for pid in probe_ids])
        padlocks = np.array([f'{gene}_{pid}' for gene in self.genes for pid in padlock_ids])
        gene_array = np.array([gene for gene in self.genes for pid in probe_ids])

        self.overlaps = generate_random_dna_sequences(12, 12, len(padlocks), True, genes=padlocks)

        #now, add overlaps to every padlock and primer, ensuring compatibility
        self.build_overlaps()

        #build the struicture for the probes dataframe
        self.probes_df = pd.DataFrame({
            'probe_id': probes,
            'gene': gene_array,
            'overlap_sequence': self.overlaps['Sequence'],
            'overlaps_reverse_complement': self.overlaps['Reverse Complement'],
            'overlaps_Tms': self.overlaps['Tms'],
            'Tm': np.zeros(len(probes)),
            'startpos': np.zeros(len(probes), dtype=int),
            'endpos': np.zeros(len(probes), dtype=int),
            'transcript_region': np.zeros(len(probes), dtype=int),
            'target': np.zeros(len(probes)),
        })
        self.probes_df.set_index('probe_id', inplace=True)
        self.probes_df['hyb'] = [gene in self.hyb_genes for gene in self.probes_df['gene']]

        # Assign targets to padlocks and primers
        self.assign_targets()

        # Assign barcodes to genes
        barcode_genes = [gene for gene in self.genes if gene not in self.hyb_genes]
        self.barcode_df = get_barcodes(self.barcode_csv_path, gene_number=len(barcode_genes))

        self.gene_dict = pd.DataFrame({'genes': barcode_genes, 'barcodes': self.barcode_df['barcode'].values, 'is_7_mer': self.barcode_df['is_7mer'].values})

        # Generate extra redundancy barcodes for all genes
        self.extra_redundancy = self.generate_extra_redundancy()

        # Assemble the final probe sequences
        self.probes_df = self.assemble_probes()

        # Save the probes dataframe, assigned barcodes, and summary
        self.save_probes(self.outdir)
        self.save_assigned_barcodes(self.outdir)
        self.write_summary(self.outdir)
        self.show_figures()

        return self.probes_df
