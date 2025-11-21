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
from ftplib import FTP, error_perm, error_temp
from pathlib import Path
from typing import Iterable, List, Mapping, Sequence
from urllib.request import urlopen

from multi_padlock_design import config

logger = logging.getLogger(__name__)

NCBI_HOST = "ftp.ncbi.nlm.nih.gov"
ENSEMBL_HOST = "ftp.ensembl.org"
GENCODE_HOST = "ftp.ebi.ac.uk"

SPECIES_CONFIG = {
    "mouse": {
        "refseq": "M_musculus/mRNA_Prot",
        "ensembl": "mus_musculus",
        "gencode": {
            "base_dir": "/pub/databases/gencode/Gencode_mouse",
            "release_prefix": "release_M",
            "filename_pattern_prefix": "gencode.vM",
        },
    },
    "human": {
        "refseq": "H_sapiens/mRNA_Prot",
        "ensembl": "homo_sapiens",
        "gencode": {
            "base_dir": "/pub/databases/gencode/Gencode_human",
            "release_prefix": "release_",
            "filename_pattern_prefix": "gencode.v",
        },
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


def _load_gencode_checksums(ftp: FTP) -> dict[str, str]:
    """Load MD5 checksums from GENCODE MD5SUMS file in current directory."""
    try:
        content = _read_remote_text_file(ftp, "MD5SUMS")
    except Exception as exc:  # pragma: no cover - network/FTP errors
        logger.warning("Unable to read GENCODE MD5SUMS file: %s", exc)
        return {}
    checksums: dict[str, str] = {}
    for line in content.splitlines():
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        checksum, filename = parts[0], parts[-1]
        checksums[filename] = checksum.lower()
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
    """Download file safely with checksum verification.

    - Writes into a temporary '.part' file per attempt.
    - Only replaces the target file after checksum passes.
    - Any error or checksum mismatch leaves the previous good file intact.
    """
    _ensure_dir(local_path.parent)

    for attempt in range(1, max_attempts + 1):
        tmp_path = local_path.with_suffix(local_path.suffix + ".part")

        # Clean up temp from previous attempts
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError as exc:
                logger.warning("Unable to remove temp file %s: %s", tmp_path, exc)

        logger.info(
            "Downloading %s to temporary file %s (attempt %d/%d)",
            filename,
            tmp_path,
            attempt,
            max_attempts,
        )

        try:
            with open(tmp_path, "wb") as handle:
                # Any network error during transfer will raise, and we will retry
                ftp.retrbinary(f"RETR {filename}", handle.write)
        except (OSError, EOFError, error_temp, error_perm) as exc:
            logger.warning(
                "Error during download of %s on attempt %d/%d: %s",
                filename,
                attempt,
                max_attempts,
                exc,
            )
            # Make sure we don't keep a partial file
            try:
                if tmp_path.exists():
                    tmp_path.unlink()
            except OSError:
                pass
            continue

        # Basic sanity check: non‑empty file
        try:
            size = tmp_path.stat().st_size
        except OSError:
            size = 0
        if size == 0:
            logger.warning(
                "Downloaded temp file %s is empty (attempt %d/%d)",
                tmp_path,
                attempt,
                max_attempts,
            )
            try:
                tmp_path.unlink()
            except OSError:
                pass
            continue

        # Checksum verification
        if _checksums_match(tmp_path, expected_checksum, checksum_type):
            # Replace final file atomically
            if local_path.exists():
                try:
                    local_path.unlink()
                except OSError as exc:
                    logger.warning(
                        "Unable to remove existing file %s before rename: %s",
                        local_path,
                        exc,
                    )
            tmp_path.replace(local_path)

            if checksum_type == "sum" and isinstance(expected_checksum, tuple):
                logger.info(
                    "%s: %d %d",
                    local_path,
                    expected_checksum[0],
                    expected_checksum[1],
                )
            logger.info("Successfully downloaded %s", local_path)
            return

        logger.warning(
            "Checksum mismatch for %s (attempt %d/%d)",
            filename,
            attempt,
            max_attempts,
        )
        # Remove broken temporary file and retry
        try:
            tmp_path.unlink()
        except OSError:
            pass

    # After all attempts we never validated checksum: fail and never keep junk
    raise RuntimeError(f"Failed to download {filename} with a valid checksum")


def _expected_extracted_path(archive_path: Path) -> Path:
    """Return the expected extracted (non-gz) path for a given archive.

    For 'foo.fa.gz' -> 'foo.fa', for 'bar.gtf.gz' -> 'bar.gtf'.
    """
    if archive_path.suffix != ".gz":
        return archive_path
    # remove only the last '.gz'
    return archive_path.with_suffix("")


def _download_files(
    ftp: FTP,
    remote_files: Iterable[str],
    destination: Path,
    *,
    force: bool = False,
    checksums: Mapping[str, object] | None = None,
    checksum_type: str | None = None,
) -> List[Path]:
    """Download archives to `destination`.

    Note: This function only knows about archives; skipping based on already
    extracted files is handled by the callers (download_refseq / ensembl /
    gencode) by filtering `remote_files` before calling this function.
    """
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
    # Archives per species go under refseq/species_key
    archives_dir = base_dir / "refseq" / species_key
    extracted_dir = base_dir / "refseq"  # extracted files live here
    _ensure_dir(archives_dir)
    if force:
        _cleanup_archives(archives_dir)

    with closing(FTP(NCBI_HOST)) as ftp:
        ftp.login()
        ftp.cwd("refseq")
        ftp.cwd(remote_subdir)
        checksum_map = _load_refseq_checksums(ftp, species_key)
        all_remote = [name for name in ftp.nlst() if "rna.fna" in name]
        if not all_remote:
            raise RuntimeError(
                f"No RefSeq files containing 'rna.fna' found in {remote_subdir}"
            )

        if not force:
            # Skip files whose *extracted* version already exists
            filtered = []
            for name in all_remote:
                archive_path = archives_dir / name
                expected = _expected_extracted_path(archive_path)
                extracted_path = extracted_dir / expected.name
                if extracted_path.exists():
                    logger.info(
                        "Skipping RefSeq %s (already have extracted %s)",
                        name,
                        extracted_path,
                    )
                    continue
                filtered.append(name)
            files_to_fetch = filtered
        else:
            files_to_fetch = all_remote

        if not files_to_fetch:
            logger.info(
                "All RefSeq files for %s already present as extracted FASTA",
                species_key,
            )
            return []

        logger.info(
            "Fetching %d RefSeq file(s) for %s", len(files_to_fetch), species_key
        )
        return _download_files(
            ftp,
            files_to_fetch,
            archives_dir,
            force=True,  # since we already filtered, force download remaining archives
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


def _latest_gencode_release(ftp: FTP, species_key: str) -> str:
    """Return the latest GENCODE release directory name for the species."""
    cfg = SPECIES_CONFIG[species_key]["gencode"]
    base_dir = cfg["base_dir"]
    release_prefix = cfg["release_prefix"]
    ftp.cwd(base_dir)
    all_entries = ftp.nlst()
    releases = [name for name in all_entries if name.startswith(release_prefix)]
    if not releases:
        raise RuntimeError(f"Unable to locate GENCODE releases under {base_dir}")

    def _release_number(name: str) -> int | None:
        """Extract numeric part from a release name.

        Examples:
        - 'release_49'  -> 49
        - 'release_M38' -> 38
        - 'release_2b'  -> None (non-numeric suffix)
        """
        tail = name.split("_", 1)[-1]  # "49", "M38", "2b", etc.
        tail = tail.lstrip("M")
        try:
            return int(tail)
        except ValueError:
            return None

    # Filter to releases with a valid numeric component
    numeric_releases = [(name, _release_number(name)) for name in releases]
    numeric_releases = [
        (name, num) for name, num in numeric_releases if num is not None
    ]

    if not numeric_releases:
        # Fall back to simple lexicographic ordering if everything is weird
        logger.warning(
            "No purely numeric GENCODE release directories found in %s; "
            "falling back to lexicographic ordering: %s",
            base_dir,
            releases,
        )
        releases.sort()
        latest = releases[-1]
    else:
        numeric_releases.sort(key=lambda x: x[1])
        latest = numeric_releases[-1][0]

    logger.info("Using GENCODE %s for %s", latest, species_key)
    return latest


def download_ensembl(
    species_key: str, base_dir: Path, *, force: bool = False
) -> tuple[str, List[Path]]:
    """Download Ensembl CDS and cDNA FASTA files for a species."""

    species_dir = SPECIES_CONFIG[species_key]["ensembl"]
    archives_base = base_dir / "ensembl" / species_key
    extracted_base = base_dir / "ensembl"
    downloaded: List[Path] = []
    if force:
        _cleanup_archives(archives_base)

    with closing(FTP(ENSEMBL_HOST)) as ftp:
        ftp.login()
        release = _latest_ensembl_release(ftp)
        for category in ("cdna", "cds"):
            remote_path = f"/pub/{release}/fasta/{species_dir}/{category}"
            ftp.cwd(remote_path)
            checksum_map = _load_ensembl_checksums(ftp)
            all_filenames = [
                name
                for name in ftp.nlst()
                if name.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz"))
            ]
            if category == "cdna":
                all_filenames = [
                    name for name in all_filenames if name.endswith(".cdna.all.fa.gz")
                ]
            if not all_filenames:
                logger.warning(
                    "No FASTA files found in %s for %s", remote_path, species_key
                )
                ftp.cwd("/")
                continue

            if not force:
                filtered: List[str] = []
                for name in all_filenames:
                    archive_path = archives_base / category / name
                    extracted_name = _expected_extracted_path(archive_path).name
                    # split long expression to satisfy line length limits
                    extracted_path = extracted_base / category / extracted_name
                    if extracted_path.exists():
                        logger.info(
                            "Skipping Ensembl %s %s (already have extracted %s)",
                            category,
                            name,
                            extracted_path,
                        )
                        continue
                    filtered.append(name)
                filenames = filtered
            else:
                filenames = all_filenames

            if not filenames:
                logger.info(
                    "All Ensembl %s files for %s already present as extracted FASTA",
                    category,
                    species_key,
                )
                ftp.cwd("/")
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
                    archives_base / category,
                    force=True,
                    checksums=checksum_map,
                    checksum_type="sum",
                )
            )
            ftp.cwd("/")
    return release, downloaded


def _http_download_with_checksum(
    url: str,
    local_path: Path,
    expected_checksum: str | None,
    max_attempts: int = 3,
) -> None:
    """Download a file over HTTP(S) with MD5 verification into `local_path`.

    Uses a temporary '.part' file and only moves it into place on success.
    """
    _ensure_dir(local_path.parent)
    for attempt in range(1, max_attempts + 1):
        tmp_path = local_path.with_suffix(local_path.suffix + ".part")
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError as exc:
                logger.warning("Unable to remove temp file %s: %s", tmp_path, exc)

        logger.info(
            "HTTP download %s -> %s (attempt %d/%d)",
            url,
            tmp_path,
            attempt,
            max_attempts,
        )
        try:
            with urlopen(url) as resp, open(tmp_path, "wb") as out:
                shutil.copyfileobj(resp, out)
        except Exception as exc:  # network errors, HTTP errors
            logger.warning(
                "Error during HTTP download of %s on attempt %d/%d: %s",
                url,
                attempt,
                max_attempts,
                exc,
            )
            try:
                if tmp_path.exists():
                    tmp_path.unlink()
            except OSError:
                pass
            continue

        try:
            size = tmp_path.stat().st_size
        except OSError:
            size = 0
        if size == 0:
            logger.warning(
                "Downloaded HTTP temp file %s is empty (attempt %d/%d)",
                tmp_path,
                attempt,
                max_attempts,
            )
            try:
                tmp_path.unlink()
            except OSError:
                pass
            continue

        if expected_checksum is None or _calculate_md5(tmp_path) == expected_checksum:
            if local_path.exists():
                try:
                    local_path.unlink()
                except OSError as exc:
                    logger.warning(
                        "Unable to remove existing file %s before rename: %s",
                        local_path,
                        exc,
                    )
            tmp_path.replace(local_path)
            logger.info("Successfully downloaded %s via HTTP", local_path)
            return

        logger.warning(
            "MD5 mismatch for HTTP download %s (attempt %d/%d)",
            url,
            attempt,
            max_attempts,
        )
        try:
            tmp_path.unlink()
        except OSError:
            pass

    raise RuntimeError(f"Failed to download {url} with a valid checksum")


def download_gencode(
    species_key: str, base_dir: Path, *, force: bool = False
) -> tuple[str, List[Path]]:
    """Download GENCODE basic annotation GTF for a species (latest release).

    Returns (release_name, [downloaded_paths]).
    """
    cfg = SPECIES_CONFIG[species_key]["gencode"]
    base_dir_cfg = cfg["base_dir"]
    filename_prefix = cfg["filename_pattern_prefix"]

    archives_dir = base_dir / "gencode" / species_key
    extracted_dir = archives_dir  # extracted GTFs live alongside per-species
    if force:
        _cleanup_archives(archives_dir)

    # First: use FTP only for discovery (release dir + filenames + MD5SUMS)
    with closing(FTP(GENCODE_HOST)) as ftp:
        ftp.login()
        release_dir = _latest_gencode_release(ftp, species_key)  # e.g. release_49
        remote_path = f"{base_dir_cfg}/{release_dir}"
        ftp.cwd(remote_path)

        checksum_map = _load_gencode_checksums(ftp)
        all_filenames = [
            name
            for name in ftp.nlst()
            if name.startswith(filename_prefix)
            and name.endswith(".basic.annotation.gtf.gz")
        ]

    if not all_filenames:
        raise RuntimeError(
            f"No GENCODE basic annotation GTF found in {remote_path} for {species_key}"
        )
    if len(all_filenames) > 1:
        logger.warning(
            "Multiple GENCODE GTF files found in %s, using the first: %s",
            remote_path,
            all_filenames[0],
        )
    chosen = all_filenames[0]

    archive_path = archives_dir / chosen
    extracted_path = extracted_dir / _expected_extracted_path(archive_path).name

    if not force and extracted_path.exists():
        logger.info(
            "Skipping GENCODE %s (already have extracted %s)",
            chosen,
            extracted_path,
        )
        return release_dir, []

    # Build HTTPS URL for the chosen file
    base_http = f"https://{GENCODE_HOST}{remote_path}"
    file_url = f"{base_http}/{chosen}"
    expected_md5 = checksum_map.get(chosen)

    logger.info(
        "Downloading GENCODE annotation via HTTPS %s for %s",
        file_url,
        species_key,
    )
    _http_download_with_checksum(file_url, archive_path, expected_md5)
    return release_dir, [archive_path]


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


def _build_gencode_rename_map(
    archives: Sequence[Path], release: str
) -> dict[Path, str]:
    """Optionally add release label to extracted GENCODE filenames."""
    if not archives:
        return {}
    # keep directory name unchanged, only adjust base file name
    rename_map: dict[Path, str] = {}
    # release is e.g. "release_49" or "release_M38"
    release_label = release
    for archive in archives:
        base_name = archive.with_suffix("").name  # drop .gz -> .gtf
        if release_label in base_name:
            new_name = base_name
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
    gencode_dir = base_dir / "gencode"

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

    gencode_release, gencode_archives = download_gencode(
        species_key, base_dir, force=args.force_download
    )
    logger.info("Downloaded %d GENCODE archive(s)", len(gencode_archives))
    gencode_rename_map = _build_gencode_rename_map(gencode_archives, gencode_release)
    gencode_extracted = _extract_archives(
        gencode_archives, gencode_dir / species_key, rename_map=gencode_rename_map
    )
    logger.info("Extracted %d GENCODE file(s)", len(gencode_extracted))


if __name__ == "__main__":
    main()
