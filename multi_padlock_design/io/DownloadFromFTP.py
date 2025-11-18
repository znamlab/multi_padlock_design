"""Download and extract reference FASTA files from NCBI RefSeq and Ensembl."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import logging
import shutil
import subprocess
from contextlib import closing
from ftplib import FTP
from pathlib import Path
from typing import Iterable, List, Mapping, Sequence

from multi_padlock_design import config

logger = logging.getLogger(__name__)

NCBI_HOST = "ftp.ncbi.nlm.nih.gov"
ENSEMBL_HOST = "ftp.ensembl.org"

SPECIES_CONFIG = {
    "mouse": {
        "refseq": "M_musculus/mRNA_Prot",
        "ensembl": "mus_musculus",
    },
    "human": {
        "refseq": "H_sapiens/mRNA_Prot",
        "ensembl": "homo_sapiens",
    },
}


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _cleanup_archives(directory: Path) -> None:
    if not directory.exists():
        return
    for gz_file in directory.rglob("*.gz"):
        try:
            gz_file.unlink()
            logger.info("Removed existing archive %s", gz_file)
        except OSError as exc:
            logger.warning("Unable to remove %s: %s", gz_file, exc)


def _read_remote_text_file(ftp: FTP, filename: str) -> str:
    buffer = io.BytesIO()
    ftp.retrbinary(f"RETR {filename}", buffer.write)
    return buffer.getvalue().decode("utf-8")


def _load_refseq_checksums(ftp: FTP, species_key: str) -> dict[str, str]:
    checksum_file = f"{species_key}.files.installed"
    try:
        content = _read_remote_text_file(ftp, checksum_file)
    except Exception as exc:  # pragma: no cover - network/FTP errors
        logger.warning("Unable to read %s checksum file: %s", checksum_file, exc)
        return {}
    checksums: dict[str, str] = {}
    for line in content.splitlines():
        parts = line.strip().split()
        if len(parts) != 2:
            continue
        checksum, filename = parts
        checksums[filename] = checksum.lower()
    return checksums


def _load_ensembl_checksums(ftp: FTP) -> dict[str, tuple[int, int]]:
    try:
        content = _read_remote_text_file(ftp, "CHECKSUMS")
    except Exception as exc:
        logger.warning("Unable to read Ensembl CHECKSUMS file: %s", exc)
        return {}
    checksums: dict[str, tuple[int, int]] = {}
    for line in content.splitlines():
        parts = line.strip().split()
        if len(parts) < 4:
            continue
        checksum_str, size_str = parts[0], parts[1]
        filename = parts[-1]
        try:
            checksum_val = int(checksum_str)
            size_val = int(size_str)
        except ValueError:
            continue
        checksums[filename] = (checksum_val, size_val)
    return checksums


def _calculate_md5(path: Path) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _calculate_sum(path: Path) -> tuple[int, int]:
    try:
        completed = subprocess.run(
            ["sum", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(
            "The 'sum' command is required to validate Ensembl downloads"
        ) from exc
    output = completed.stdout.strip().split()
    if len(output) < 2:
        raise RuntimeError(f"Unexpected sum output for {path}: {completed.stdout}")
    return int(output[0]), int(output[1])


def _checksums_match(
    path: Path,
    expected: object | None,
    checksum_type: str | None,
) -> bool:
    if not expected or not checksum_type:
        return True
    if checksum_type == "md5":
        return _calculate_md5(path) == str(expected).lower()
    if checksum_type == "sum":
        if not isinstance(expected, tuple) or len(expected) != 2:
            return False
        actual = _calculate_sum(path)
        return actual == expected
    raise ValueError(f"Unsupported checksum type '{checksum_type}'")


def _download_with_retries(
    ftp: FTP,
    filename: str,
    local_path: Path,
    expected_checksum: object | None,
    checksum_type: str | None,
    max_attempts: int = 3,
) -> None:
    for attempt in range(1, max_attempts + 1):
        with open(local_path, "wb") as handle:
            ftp.retrbinary(f"RETR {filename}", handle.write)
        if _checksums_match(local_path, expected_checksum, checksum_type):
            if checksum_type == "sum" and isinstance(expected_checksum, tuple):
                logger.info(
                    "%s: %d %d",
                    local_path,
                    expected_checksum[0],
                    expected_checksum[1],
                )
            return
        logger.warning(
            "Checksum mismatch for %s (attempt %d/%d)",
            local_path,
            attempt,
            max_attempts,
        )
    raise RuntimeError(f"Failed to download {filename} with a valid checksum")


def _download_files(
    ftp: FTP,
    remote_files: Iterable[str],
    destination: Path,
    *,
    force: bool = False,
    checksums: Mapping[str, object] | None = None,
    checksum_type: str | None = None,
) -> List[Path]:
    downloaded: List[Path] = []
    _ensure_dir(destination)
    for filename in remote_files:
        local_path = destination / filename
        expected = checksums.get(filename) if checksums else None
        needs_download = force or not local_path.exists()
        if not needs_download and expected:
            if not _checksums_match(local_path, expected, checksum_type):
                logger.warning(
                    "Existing file %s has invalid checksum; redownloading",
                    local_path,
                )
                needs_download = True
        if needs_download:
            logger.info("Fetching %s", filename)
            _download_with_retries(
                ftp,
                filename,
                local_path,
                expected,
                checksum_type,
            )
            logger.info("Stored %s", local_path)
        elif expected:
            if checksum_type == "sum" and isinstance(expected, tuple):
                logger.info(
                    "%s: %d %d (verified)",
                    local_path,
                    expected[0],
                    expected[1],
                )
            else:
                logger.info("Verified existing file %s", local_path)
        else:
            logger.info("Using existing file %s", local_path)
        downloaded.append(local_path)
    return downloaded


def download_refseq(
    species_key: str, base_dir: Path, *, force: bool = False
) -> List[Path]:
    """Download RefSeq transcript FASTA files for a species."""

    remote_subdir = SPECIES_CONFIG[species_key]["refseq"]
    output_dir = base_dir / "refseq" / species_key
    _ensure_dir(output_dir)
    if force:
        _cleanup_archives(output_dir)

    with closing(FTP(NCBI_HOST)) as ftp:
        ftp.login()
        ftp.cwd("refseq")
        ftp.cwd(remote_subdir)
        checksum_map = _load_refseq_checksums(ftp, species_key)
        files_to_fetch = [name for name in ftp.nlst() if "rna.fna" in name]
        if not files_to_fetch:
            raise RuntimeError(
                f"No RefSeq files containing 'rna.fna' found in {remote_subdir}"
            )
        logger.info(
            "Fetching %d RefSeq file(s) for %s", len(files_to_fetch), species_key
        )
        return _download_files(
            ftp,
            files_to_fetch,
            output_dir,
            force=force,
            checksums=checksum_map,
            checksum_type="md5",
        )


def _latest_ensembl_release(ftp: FTP) -> str:
    ftp.cwd("/pub")
    releases = [name for name in ftp.nlst() if name.startswith("release-")]
    if not releases:
        raise RuntimeError("Unable to locate Ensembl release directories")
    releases.sort(key=lambda value: int(value.split("-", 1)[1]))
    latest = releases[-1]
    logger.info("Using Ensembl %s", latest)
    return latest


def download_ensembl(
    species_key: str, base_dir: Path, *, force: bool = False
) -> tuple[str, List[Path]]:
    """Download Ensembl CDS and cDNA FASTA files for a species."""

    species_dir = SPECIES_CONFIG[species_key]["ensembl"]
    output_dir = base_dir / "ensembl" / species_key
    downloaded: List[Path] = []
    if force:
        _cleanup_archives(output_dir)

    with closing(FTP(ENSEMBL_HOST)) as ftp:
        ftp.login()
        release = _latest_ensembl_release(ftp)
        for category in ("cdna", "cds"):
            remote_path = f"/pub/{release}/fasta/{species_dir}/{category}"
            ftp.cwd(remote_path)
            checksum_map = _load_ensembl_checksums(ftp)
            filenames = [
                name
                for name in ftp.nlst()
                if name.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz"))
            ]
            if category == "cdna":
                filenames = [
                    name for name in filenames if name.endswith(".cdna.all.fa.gz")
                ]
            if not filenames:
                logger.warning(
                    "No FASTA files found in %s for %s", remote_path, species_key
                )
                continue
            logger.info(
                "Fetching %d Ensembl %s file(s) for %s",
                len(filenames),
                category,
                species_key,
            )
            downloaded.extend(
                _download_files(
                    ftp,
                    filenames,
                    output_dir / category,
                    force=force,
                    checksums=checksum_map,
                    checksum_type="sum",
                )
            )
            ftp.cwd("/")
    return release, downloaded


def _extract_archives(
    archives: Sequence[Path],
    destination: Path,
    *,
    rename_map: Mapping[Path, str] | None = None,
) -> List[Path]:
    extracted: List[Path] = []
    if not archives:
        return extracted
    _ensure_dir(destination)
    for archive in archives:
        if archive.suffix != ".gz":
            logger.debug("Skipping non-gzip file %s", archive)
            continue
        target_name = (
            rename_map[archive]
            if rename_map and archive in rename_map
            else archive.with_suffix("").name
        )
        target_path = destination / target_name
        if target_path.exists():
            logger.info("Replacing existing file %s", target_path)
            target_path.unlink()
        with gzip.open(archive, "rb") as src, open(target_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        archive.unlink()
        extracted.append(target_path)
        logger.info("Extracted %s", target_path)
    return extracted


def _normalize_release_label(release: str) -> str:
    return release if release.startswith("release-") else f"release-{release}"


def _build_ensembl_rename_map(
    archives: Sequence[Path], release: str
) -> dict[Path, str]:
    if not archives:
        return {}
    release_label = _normalize_release_label(release)
    rename_map: dict[Path, str] = {}
    for archive in archives:
        base_name = archive.with_suffix("").name
        if ".cdna." in base_name:
            new_name = base_name.replace(".cdna.", f".{release_label}.cdna.", 1)
        elif ".cds." in base_name:
            new_name = base_name.replace(".cds.", f".{release_label}.cds.", 1)
        else:
            new_name = f"{base_name}.{release_label}"
        rename_map[archive] = new_name
    return rename_map


def _resolve_species(species_override: str | None = None) -> str:
    species_value = (species_override or getattr(config, "species", "")).strip().lower()
    if not species_value:
        raise RuntimeError("config.species must be set to 'mouse' or 'human'")
    if species_value not in SPECIES_CONFIG:
        raise ValueError(
            f"Unsupported species '{species_override or config.species}'. "
            f"Supported: {list(SPECIES_CONFIG)}"
        )
    return species_value


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download RefSeq and Ensembl transcript FASTA files and place the "
            "extracted FASTA files under config.BASE_DIR/refseq and /ensembl. "
            "Species defaults to config.species (mouse or human)."
        )
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Redownload and re-extract all archives even if files already exist.",
    )
    parser.add_argument(
        "--species",
        choices=sorted(SPECIES_CONFIG.keys()),
        help=(
            "Override config.species for this run. Choices: %(choices)s. "
            "Defaults to the species declared in config.species."
        ),
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = _parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    species_key = _resolve_species(args.species)
    base_dir = Path(config.BASE_DIR)
    logger.info("Saving files under %s", base_dir)

    refseq_dir = base_dir / "refseq"
    ensembl_dir = base_dir / "ensembl"

    refseq_archives = download_refseq(species_key, base_dir, force=args.force_download)
    logger.info("Downloaded %d RefSeq archive(s)", len(refseq_archives))
    refseq_extracted = _extract_archives(refseq_archives, refseq_dir)
    logger.info("Extracted %d RefSeq file(s)", len(refseq_extracted))

    ensembl_release, ensembl_archives = download_ensembl(
        species_key, base_dir, force=args.force_download
    )
    logger.info("Downloaded %d Ensembl archive(s)", len(ensembl_archives))
    ensembl_rename_map = _build_ensembl_rename_map(ensembl_archives, ensembl_release)
    ensembl_extracted = _extract_archives(
        ensembl_archives, ensembl_dir, rename_map=ensembl_rename_map
    )
    logger.info("Extracted %d Ensembl file(s)", len(ensembl_extracted))


if __name__ == "__main__":
    main()
