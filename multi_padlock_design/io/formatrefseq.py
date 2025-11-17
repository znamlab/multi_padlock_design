# fasta_db.py

import gzip
import os
import re
import subprocess
from typing import Iterable, List, Tuple, Union

import multi_padlock_design.config as config

# ---------------------------
# FASTA reading (gz or plain)
# ---------------------------


def _open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def _iter_fasta(paths: Iterable[str]):
    """Yield (header, seq) from one or many FASTA files (handles multiline seqs)."""
    header = None
    seq_chunks = []
    for path in paths:
        with _open_text(path) as f:
            for line in f:
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(seq_chunks)
                    header = line.rstrip("\n")
                    seq_chunks = []
                else:
                    seq_chunks.append(line.strip())
    if header is not None:
        yield header, "".join(seq_chunks)


# ---------------------------
# Header detection & parsing
# ---------------------------

_re_paren = re.compile(r"\(([^)]+)\)")
_re_gene_symbol = re.compile(r"gene_symbol:([^\s]+)")


def _detect_source(header: str) -> str:
    h = header[1:] if header.startswith(">") else header
    if h.startswith(("NM_", "NR_")) or "|NM_" in h or "|NR_" in h:
        return "refseq"
    if "gene_symbol:" in h or h.startswith(("ENST", "ENSMUST", "ENS")):
        return "ensembl"
    if "|" in h:
        return "custom"
    return "unknown"


def _acronym_from_header(header: str, source: str = "auto") -> str:
    """Return a gene symbol/acronym from various header styles."""
    h = header[1:] if header.startswith(">") else header
    src = _detect_source(header) if source == "auto" else source

    if src == "refseq":
        # e.g. >NM_... Mus musculus FMR1 ... (Fxr1), transcript ...
        paren = _re_paren.findall(h)
        if paren:
            return paren[-1]
        # fallback to first token after accession
        toks = h.split()
        return toks[1] if len(toks) > 1 else toks[0]

    if src == "ensembl":
        # e.g. ... gene_symbol:Trav3-1 ...
        m = _re_gene_symbol.search(h)
        if m:
            return m.group(1)
        # fallback: last (...) group if any
        paren = _re_paren.findall(h)
        if paren:
            return paren[-1]
        # last resort: first token
        return h.split()[0]

    if src == "custom":
        # e.g. >Olfr1366|OTTMUSG... -> take the first pipe-delimited token
        return h.split("|", 1)[0]

    # Unknown: try (), else first token
    paren = _re_paren.findall(h)
    if paren:
        return paren[-1]
    return h.split()[0]


def _keep_record(header: str, source: str, keep_nm_nr_only: bool) -> bool:
    """Filter logic (preserves your old behavior for RefSeq)."""
    src = _detect_source(header) if source == "auto" else source
    if src == "refseq" and keep_nm_nr_only:
        acc = header.split()[0].lstrip(">")
        return acc.startswith(("NM_", "NR_"))
    # For Ensembl/custom we keep all by default
    return True


# ---------------------------
# Main routine
# ---------------------------


def fastadb(
    indir: str,
    files: Union[Tuple[str, str, int], List[str]],
    name: str,
    source: str = "auto",
    keep_nm_nr_only: bool = True,
) -> Tuple[List[str], List[str], List[str], str]:
    """
    Build helper files and a consolidated selected FASTA.

    Args:
        indir: directory containing the FASTA files
        files:
          - (prefix, suffix, filenum)  -> builds prefix{1..N}{suffix}
          - or a list of filenames     -> uses exactly these
        name: base name for outputs
        source: 'auto'|'refseq'|'ensembl'|'custom'
        keep_nm_nr_only: only applies to RefSeq; keeps NM/NR like the original code

    Returns:
        (all_headers, selected_headers, acronyms, selected_fasta_path)
    """
    if isinstance(files, tuple):
        prefix, suffix, n = files
        paths = [os.path.join(indir, f"{prefix}{i+1}{suffix}") for i in range(n)]
    else:
        paths = [os.path.join(indir, p) for p in files]

    # Collect, and simultaneously write “all” files
    all_headers: List[str] = []
    all_seqs: List[str] = []
    sel_headers: List[str] = []
    sel_seqs: List[str] = []
    acronyms: List[str] = []

    for hdr, seq in _iter_fasta(paths):
        all_headers.append(hdr)
        all_seqs.append(seq)
        if _keep_record(hdr, source, keep_nm_nr_only):
            sel_headers.append(hdr)
            sel_seqs.append(seq)
            acronyms.append(_acronym_from_header(hdr, source))

    # Write outputs
    allh_path = os.path.join(indir, f"{name}.allheaders.txt")
    alls_path = os.path.join(indir, f"{name}.allseqs.txt")
    selh_path = os.path.join(indir, f"{name}.selectedheaders.txt")
    sels_path = os.path.join(indir, f"{name}.selectedseqs.txt")
    acr_path = os.path.join(indir, f"{name}.acronymheaders.txt")
    selfa_path = os.path.join(indir, f"{name}.selected.fasta")

    with open(allh_path, "w") as fh, open(alls_path, "w") as fs:
        for h, s in zip(all_headers, all_seqs):
            fh.write(h + "\n")
            fs.write(s + "\n")

    with (
        open(selh_path, "w") as fh,
        open(sels_path, "w") as fs,
        open(selfa_path, "w") as ff,
    ):
        for h, s in zip(sel_headers, sel_seqs):
            fh.write(h + "\n")
            fs.write(s + "\n")
            # also write a consolidated FASTA that BLAST can read directly
            ff.write(h + "\n")
            # wrap to 60 chars for nicer FASTA (optional)
            for i in range(0, len(s), 60):
                ff.write(s[i : i + 60] + "\n")

    with open(acr_path, "w") as fa:
        for a in acronyms:
            fa.write(a + "\n")

    return all_headers, sel_headers, acronyms, selfa_path


# ---------------------------
# BLAST helpers
# ---------------------------


def blastdb(species: str = "mouse", dbtype: str = "nucl"):
    """Create a BLAST database from a single FASTA file."""
    fasta_path = config.blast_db_file
    if not os.path.isfile(fasta_path):
        print(f"FASTA file not found: {fasta_path}, attempting to create one.")
        # Build the selected FASTA directly, avoiding retrieveseq.loaddb to break the cycle
        fastadir = os.path.join(config.BASE_DIR, config.reference_transcriptome)
        try:
            if getattr(config, "reference_transcriptome", "").lower() == "refseq":
                pre_suf = getattr(
                    config, "fasta_pre_suffix_refseq", None
                )  # expected (prefix, suffix)
                filenum = getattr(config, "fasta_filenum_refseq", None)  # expected int
                if (
                    isinstance(pre_suf, (list, tuple))
                    and len(pre_suf) == 2
                    and isinstance(filenum, int)
                ):
                    fastadb(
                        fastadir,
                        (pre_suf[0], pre_suf[1], filenum),
                        species,
                        source="refseq",
                        keep_nm_nr_only=True,
                    )
                else:
                    raise ValueError(
                        "Invalid RefSeq config: expected fasta_pre_suffix_refseq=(prefix, suffix) and fasta_filenum_refseq=int"
                    )
            elif getattr(config, "reference_transcriptome", "").lower() == "ensembl":
                files = [config.cdna_file]
                if hasattr(config, "extra_files") and config.extra_files:
                    files += config.extra_files
                fastadb(
                    fastadir,
                    files,
                    species,
                    source="ensembl",
                    keep_nm_nr_only=False,
                )
            else:
                # Fallback: try to build with whatever cdna_file points to
                files = [getattr(config, "cdna_file", fasta_path)]
                fastadb(fastadir, files, species)
        except Exception as e:
            print(f"Failed to prepare FASTA for BLAST: {e}")
            return

    try:
        # Use the transcriptome directory for BLAST outputs
        transcriptome_dir = os.path.join(
            config.BASE_DIR, config.reference_transcriptome
        )
        out_prefix = os.path.join(transcriptome_dir, species + ".transcriptome")
        # Avoid re-creating if DB already present (check a couple of index files)
        if not (
            os.path.isfile(out_prefix + ".nal") or os.path.isfile(out_prefix + ".nhr")
        ):
            subprocess.run(
                [
                    "makeblastdb",
                    "-in",
                    fasta_path,
                    "-dbtype",
                    dbtype,
                    "-out",
                    out_prefix,
                    # "-title", f'"{species}_transcriptome"',
                ],
                check=True,
            )
        else:
            print(f"BLAST database already exists for {species} at {out_prefix}.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while creating BLAST database: {e}")
