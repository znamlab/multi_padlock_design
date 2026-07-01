# Probe construction (`makesnails.py`)

This document describes how `snailmaker` (in [`lib/makesnails.py`](../lib/makesnails.py))
turns a set of BLAST-verified candidate targets into a full panel of SNAIL
(padlock + primer) probes, ready to order.

## 1. Background: what a SNAIL probe pair is

Each gene in the panel is detected by 1-4 independent **padlock/primer pairs**
("SNAIL" probes). Padlocks and primers work together:

- The **primer** binds the RNA/cDNA target immediately downstream of the padlock's
  binding site and, once extended, provides the substrate the padlock ligates onto.
- The **padlock** binds the target on both sides of a nick, is circularized by
  ligation, and is then amplified by rolling circle amplification (RCA) to produce
  a bright, localized signal.
- Padlock and primer share a common **overlap sequence**, split across their two
  ends, so that only a matched padlock/primer pair for the *same* target site
  ligate/prime correctly and cross-reactivity between unrelated probes is
  suppressed.
- Downstream, the RCA product is read out either:
  - **by hybridization** — a single fixed fluorescent probe (a "hyb probe")
    binds a fixed sequence built into the padlock (used for a small set of
    genes with dedicated imaging channels), or
  - **by sequencing** — a gene-specific **barcode** plus a per-pair **UMI**
    built into the padlock is read out over several sequencing-by-ligation/
    hybridization cycles (the standard case for most genes in the panel).

`snailmaker` automates picking targets, generating the overlap sequences,
assigning barcodes, and assembling the final oligo sequences for both readout
schemes.

## 2. Inputs

| Input | Purpose | Loaded by |
|---|---|---|
| `outdir/panel_name/<gene>/3.AllSpecificTargets_*.csv` | Per-gene BLAST-verified candidate targets (one subdirectory per gene) | `read_candidates` |
| `hyb_probes_path` (csv) | Available hybridization readout probes, columns include `SNAIL TARG` (gene) and `SEQ` (probe sequence, `/`-delimited fields) | loaded into `self.hyb_probes` |
| `barcode_csv_path` (csv) | Pool of available sequencing barcodes, columns `barcode`, `is_7mer`, `gene`, `panel` (see `get_barcodes`) | `get_barcodes` |
| `hyb_genes` (list) | Genes to read out by hybridization instead of sequencing | constructor arg |
| `min_probes` (int, ≤4) | Padlock/primer pairs per gene; genes with fewer usable candidates are dropped | constructor arg |
| `sequencing_primer` (str) | Constant sequencing-primer sequence built into every non-hyb padlock | constructor arg (default `GTCGGACTGTAGAACTCT`) |

`min_probes` is capped at 4 because each pair is tagged with a single-base
UMI (`A`/`T`/`C`/`G`) — there's no fifth distinguishable UMI value.

## 3. Candidate target selection

- `read_candidates(outdir, panel_name, probefile)` walks every gene
  subdirectory under `outdir/panel_name`, reads its
  `3.AllSpecificTargets_*.csv` (a semicolon-safe custom csv — commas inside
  `[...]` are protected by `smart_split`), tags each row with its gene name,
  and concatenates everything into `self.big_df`. Rows with no `acronym`
  are dropped. This runs automatically in `__init__`.
- `probes_per_gene()` plots the number of candidate targets per gene
  (sorted descending) with a red line at `min_probes`, so you can see at a
  glance which genes are at risk of being dropped. The figure is stashed in
  `self.figures` and shown later by `show_figures()` if `plot_report=True`.
- `choose_targets(gene, n_probes)` filters `big_df` to one gene and ranks
  candidates by a **homopolymer score** (2 = no run of 4+ identical bases,
  1 = has a run of 4, 0 = has a run of 5+), preferring probes least likely to
  suffer synthesis/sequencing errors from long homopolymers. It returns the
  top `n_probes` rows.
- `assign_targets()` calls `choose_targets` for every gene in `self.genes`
  and writes each chosen target's `target`/`Tm`/`startpos`/`endpos`/
  `transcript_region` into both the padlock and primer row for that pair in
  `self.probes_df`.
  - A gene with **zero** candidates, or **fewer candidates than
    `min_probes`**, is dropped entirely: its rows are removed from
    `self.probes_df`, its name removed from `self.genes`, and it's recorded
    in `self.dropped_genes`. The console warning and the final
    `summary.txt` both flag this in capitals — `... gene DROPPED`.

## 4. Overlap sequence generation

- `generate_random_dna_sequences(min_barcode, length, num_sequences, for_overlap, genes)`
  is the general-purpose random-oligo generator used for the padlock/primer
  overlap. It keeps drawing random sequences until it has `num_sequences`
  that all satisfy:
  - melting temperature (nearest-neighbour) between 25-35°C,
  - no homopolymer run of 3+,
  - (if `for_overlap`) bases `[3:6]` don't produce a ligation-blocking motif
    at the junction, and
  - at least 2 mismatches against every other sequence's first
    `min_barcode` bases, and at least 4 mismatches against every other
    sequence overall — so overlaps stay distinguishable even under
    sequencing error.
- `main()` calls this once per **padlock** (12bp, `min_barcode=12`,
  `for_overlap=True`), producing one unique overlap sequence per padlock.
- `build_overlaps()` then duplicates each padlock's overlap row into a
  matching **primer** row (same sequence, ID offset by `+10`, e.g.
  `GENE_1` → `GENE_11`), so every padlock/primer pair shares one overlap
  sequence split into a 5' half and 3' half (`overlap_5_seq` /
  `overlap_3_seq`, see §6).

## 5. Barcode & redundancy assignment

- Genes in `hyb_genes` are read out by hybridization and don't need a
  sequencing barcode. All other genes (`barcode_genes`) get one:
  `get_barcodes(barcode_csv_path, gene_number=len(barcode_genes))` pulls the
  first N **unused** barcodes (`gene` column blank) from the barcode pool,
  in file order — the barcode CSV is expected to already be sorted so that
  preferred barcodes (e.g. 7-mers) come first. This build's result
  (`self.gene_dict`, mapping `genes` → `barcodes`/`is_7_mer`) is written out
  per-panel as `assigned_barcodes.csv` (see §7) — **the source barcode CSV
  itself is never modified**, so barcode usage must be tracked across panels
  by inspecting each panel's `assigned_barcodes.csv`.
- `generate_extra_redundancy()` gives every gene (hyb and non-hyb alike) a
  random homopolymer-free 5-mer, used as a short extra identifying/padding
  sequence in the final padlock (see §6). This is independent of the
  sequencing barcode.

## 6. Final sequence assembly (`assemble_probes`)

For every gene, for every pair (`1..min_probes`), the padlock and primer
sequences are built 5'→3' as:

**Primer** (same for hyb and sequencing readout):

```
target_seq + overlap_sequence
```

**Padlock, hybridization readout** (`gene in hyb_genes`):

```
overlap_5' + target_seq + hyb_sequence + extra_redundancy[0:2] + overlap_3'
```

**Padlock, sequencing readout** (all other genes):

```
overlap_5' + target_seq + sequencing_primer + barcode + UMI + extra_redundancy + overlap_3'
```

Where:
- `target_seq` is the reverse complement of the first 22 bases of the chosen
  target (padlock) or of bases `[23:]` (primer) — the 1-base gap between the
  two halves is intentional, leaving room at the ligation junction.
- `overlap_5'` / `overlap_3'` are the reverse complements of the first 6 and
  remaining bases of the shared 12bp overlap sequence (§4).
- `hyb_sequence` is looked up from `hyb_probes` by gene name, reverse
  complemented.
- `barcode` and `is_7_mer` come from `self.gene_dict`; `UMI` is a single base
  (`A`/`T`/`C`/`G`) selected by pair index, so each of a gene's up-to-4 pairs
  gets a distinct UMI.
- `extra_redundancy` is the gene's 5-mer from §5 (only its first 2 bases are
  used in the hyb case, to keep the RCA product length in the target range).

The assembled sequences and their lengths are written back into
`self.probes_df` as `full_sequence` / `length`.

## 7. Outputs

Calling `snailmaker.main()` runs the whole pipeline and writes three files to
`outdir/panel_name/`:

| File | Contents |
|---|---|
| `probes_df.csv` | One row per padlock/primer, including the final `full_sequence` |
| `assigned_barcodes.csv` | `self.gene_dict` — the barcode/UMI-pool assignment made for this panel run |
| `summary.txt` | Gene count & list, hyb probe count, hyb-readout genes, barcode count, 7-mer count, and the list of `DROPPED GENES` |

If `plot_report=True`, `show_figures()` is also called at the end, displaying
every figure collected during the run (currently just `probes_per_gene`).

## 8. Pipeline order (what `main()` does)

```
1. probes_per_gene()               -- QC plot of candidates/gene
2. build padlock/primer overlap sequences (generate_random_dna_sequences)
3. build_overlaps()                -- mirror overlaps onto primer rows
4. build empty probes_df skeleton  -- one row per padlock/primer
5. assign_targets()                -- pick + assign targets, drop under-supported genes
6. assign barcodes to non-hyb genes -> self.gene_dict
7. generate_extra_redundancy()     -- per-gene 5-mer
8. assemble_probes()               -- build full_sequence for every probe
9. save_probes / save_assigned_barcodes / write_summary
10. show_figures()
```

Each step can also be called individually for more control (e.g. to inspect
`self.big_df` or `self.probes_df` between steps, or to re-run just
`assemble_probes()` after manually editing `self.gene_dict`).
