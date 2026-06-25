# Input file formats

The pipeline accepts two kinds of input (chosen via the `csv` / `fasta` prompt or
the `INPUT_TYPE` SLURM variable):

1. **CSV gene list** — gene symbols, optionally carrying per-gene barcode and
   linker (backbone) sequences. Described in detail below.
2. **FASTA file** — one target sequence per record. See [FASTA input](#fasta-input).

---

## CSV gene list

A plain comma-separated text file. **One gene per row, no header line.**

The columns are **positional** (there are no column names). Only the first
column is required; the remaining columns supply the padlock backbone and are
optional.

| Position | Column              | Required | Default if omitted/empty           | Used as            |
|---------:|---------------------|:--------:|------------------------------------|--------------------|
| 1        | `gene`              | yes      | —                                  | gene symbol/acronym lookup |
| 2        | `linker_upstream`   | no       | `ATCGTCGGACTGTAGAACTCTGAACCTGTCG`  | `linker1` (backbone, 5′ of barcode) |
| 3        | `barcode`           | no       | `NNNNNNNNNN`                       | barcode            |
| 4        | `linker_downstream` | no       | `CT`                               | `linker2` (backbone, 3′ of barcode) |

Any columns beyond the fourth are ignored.

Parsed by [`readgenefile`](../multi_padlock_design/io/checkinput.py); columns 2–4
are consumed by
[`assembleprobes`](../multi_padlock_design/select_probes/finalizeprobes.py).

### How the columns become a padlock

For each accepted target site, a `2 × arm_length` target-complementary sequence
(`probe`) is split in half and the backbone is inserted between the halves. The
assembled padlock, written 5′→3′, is:

```
5'Phos − [arm B] − linker_upstream − barcode − linker_downstream − [arm A] − 3'
            │            │              │              │              │
   probe[arm_len:]   column 2       column 3       column 4    probe[0:arm_len]
```

- `arm B` = `probe[arm_length:]` (carries the ligatable 5′-phosphate end)
- `arm A` = `probe[0:arm_length]` (carries the 3′-OH end)
- `linker_upstream`, `barcode`, `linker_downstream` form the connecting loop, in
  exactly that order — hence the column names.

The `/5Phos/` modification is **not** added by the code; it is a synthesis
modification specified when ordering, and is functionally required because
circularization ligates arm B's 5′-phosphate to arm A's 3′-OH.

### Rules and gotchas

- **Positional, not named.** To set `linker_downstream` you must also fill
  columns 2 and 3, otherwise your downstream linker would be read as the
  barcode. Use the default values explicitly if you only want to override a
  later column.
- **Empty == default.** A missing column *or* an empty field (e.g. a trailing
  comma) falls back to the default in the table above. You cannot produce a
  zero-length linker/barcode this way.
- **No spaces around commas.** The gene list is parsed with a plain
  `line.split(",")` (not a quoting CSV parser), and fields are **not** trimmed.
  `Gad1, ACCC` would make the barcode `" ACCC"` (with a leading space). Write
  fields back-to-back: `Gad1,...,ACCC,...`.
- **Sequences are concatenated verbatim.** Provide DNA bases (`ACGT`) for
  linker/barcode columns in the exact case you want; they are not validated,
  uppercased, or screened for secondary structure during assembly.
- **Duplicate gene names are dropped**, keeping the first occurrence (and its
  linkers).
- **Bare gene-only rows are valid.** A file with just gene symbols gives every
  probe the default backbone and an `NNNNNNNNNN` placeholder barcode — fine when
  you are not barcoding, but not a usable barcode set.

### Examples

**Bare gene list** (default backbone, placeholder barcode for every gene):

```
Gad1
Gad2
Slc17a7
```

**Full backbone control** (all four columns; this is the layout that produces
the probe structure shown above):

```
Gad1,ATCGTCGGACTGTAGAACTCTGAACCTGTCG,AGGGCATGTA,CCT
Gad2,ATCGTCGGACTGTAGAACTCTGAACCTGTCG,CTCAACGAAG,CCT
Slc17a7,ATCGTCGGACTGTAGAACTCTGAACCTGTCG,GTGCACATAG,CCT
```

Here each gene keeps the default upstream linker but uses a unique 10-nt barcode
and a custom 3-nt downstream linker (`CCT`).

See also [examples/input_human_genes.txt](../examples/input_human_genes.txt)
(bare) and
[examples/input_human_genes_withbackbone.csv](../examples/input_human_genes_withbackbone.csv)
(with backbone).

### SLURM note

`slurm-probe-design` first splits the CSV into one single-row file per gene via
[`csv_splitter`](../multi_padlock_design/io/csv_splitter.py), naming each file
`<gene>.csv`. All columns are preserved in the split files, so per-gene barcodes
and linkers carry through to the per-gene jobs.

---

## FASTA input

A FASTA file with one target sequence per record. Each header becomes the gene
name (sanitized of characters such as `:/\[]?"<> ` and commas). When FASTA input
is used, every record receives the **default** backbone and barcode — there is
no per-record barcode/linker mechanism in this mode.

---

## Per-probe identifiers (optional)

By default every padlock for a gene shares the same backbone, so individual
probes of a gene are indistinguishable. Enabling per-probe identifiers inserts a
short, error-robust code **between the barcode and the downstream linker**, so
each padlock within a gene is uniquely labelled:

```
… barcode − IDENTIFIER − linker_downstream − [arm A] − 3'
```

Enable it in [config.py](../multi_padlock_design/config.py):

- `use_padlock_identifiers = True`
- `padlock_identifier_file` — defaults to
  [barcodes/padlock_identifiers_3bp.txt](../multi_padlock_design/barcodes/padlock_identifiers_3bp.txt)
- `padlock_identifier_min_distance` — minimum Hamming distance of the robust tier
  (default 2)

Behaviour:

- Applied to the **final subset only** (`6.ProbesRandomSubsetN=*`), which gains a
  `padlock_id` column recording each probe's code. Requires a per-gene probe
  count ("Number of probes per gene").
- The identifier file is an **ordered, tiered** list. The leading *robust tier*
  is pairwise ≥ `min_distance` apart (the 3bp set: 16 codes, single-error
  detecting) and is ordered so the first base distinguishes the first 4 probes,
  the first two bases distinguish all 16, and the third base is the error check.
- Beyond the robust tier the list is **extended** with distinct (distance-1)
  codes (the 3bp set reaches 60 total). A gene that uses these emits a warning
  (still unique, but not error-robust); a gene that needs more codes than the
  file provides raises an error.

Generate other sets with
[`padlock_identifiers.py`](../multi_padlock_design/barcodes/padlock_identifiers.py):

```
python -m multi_padlock_design.barcodes.padlock_identifiers \
    --length 4 --min-distance 3 --out my_ids_4bp.txt
```
