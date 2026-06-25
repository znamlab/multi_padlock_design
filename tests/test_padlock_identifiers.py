"""Tests for per-probe (padlock-specific) identifiers."""

import pytest

from multi_padlock_design.barcodes import padlock_identifiers as pi
from multi_padlock_design.io import createoutput
from multi_padlock_design.select_probes import finalizeprobes


def _hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


# --------------------------------------------------------------------------- #
# Identifier set: ordering, tiering, validation
# --------------------------------------------------------------------------- #
def test_committed_3bp_file_properties():
    codes = pi.load_identifiers(_file_path())
    # tiered robust prefix
    assert pi.robust_prefix_size(codes, 2) == 16
    # informative ordering
    assert len({c[0] for c in codes[:4]}) == 4  # first base separates first 4
    assert len({c[:2] for c in codes[:16]}) == 16  # first 2 bases separate first 16
    # robust tier is pairwise >= distance 2
    for i in range(16):
        for j in range(i + 1, 16):
            assert _hamming(codes[i], codes[j]) >= 2
    # whole list distinct and homopolymer-free
    assert len(set(codes)) == len(codes)
    assert all(pi.homopolymer_free(c) for c in codes)


def test_informative_order_cycles_first_base():
    alphabetical = sorted(
        [
            "AAG",
            "ACT",
            "AGC",
            "ATA",
            "CAC",
            "CCA",
            "CGG",
            "CTT",
            "GAA",
            "GCC",
            "GGT",
            "GTG",
            "TAT",
            "TCG",
            "TGA",
            "TTC",
        ]
    )
    ordered = pi.informative_order(alphabetical)
    assert len({c[0] for c in ordered[:4]}) == 4
    assert len({c[:2] for c in ordered[:16]}) == 16
    assert sorted(ordered) == alphabetical  # same set, just reordered


def test_extend_distinct_keeps_robust_prefix_and_dedupes():
    robust = ["AAG", "CCA"]
    pool = pi.candidate_pool(3, exclude_homopolymers=True)
    extended = pi.extend_distinct(robust, pool)
    assert extended[: len(robust)] == robust  # robust prefix preserved
    assert len(set(extended)) == len(extended)  # all distinct
    assert set(robust).issubset(extended)


def test_robust_prefix_size():
    assert pi.robust_prefix_size(["AAA", "CCC", "AAC"], 2) == 2  # AAC vs AAA == 1
    assert pi.robust_prefix_size(["AAG", "CCA", "GTG"], 2) == 3


def test_validate_rejects_duplicates_and_unequal_length():
    with pytest.raises(ValueError):
        pi.validate_identifier_set(["AAG", "AAG"], 2)
    with pytest.raises(ValueError):
        pi.validate_identifier_set(["AAG", "AC"], 2)
    with pytest.raises(ValueError):
        pi.validate_identifier_set([], 2)


# --------------------------------------------------------------------------- #
# Assembly
# --------------------------------------------------------------------------- #
def _genepars(linkers):
    # assembleprobes only uses genepars[1]
    return (None, linkers, None)


def test_assembleprobes_inserts_identifier_between_barcode_and_linker():
    targets = [["ACGTTGCA"]]  # armlength 4 -> armB="TGCA", armA="ACGT"
    genepars = _genepars([["LU", "BC", "LD"]])
    out = finalizeprobes.assembleprobes(targets, genepars, 4, identifiers=[["ID"]])
    assert out == [["TGCA" + "LU" + "BC" + "ID" + "LD" + "ACGT"]]


def test_assembleprobes_unchanged_without_identifiers():
    targets = [["ACGTTGCA"]]
    genepars = _genepars([["LU", "BC", "LD"]])
    out = finalizeprobes.assembleprobes(targets, genepars, 4)
    assert out == [["TGCA" + "LU" + "BC" + "LD" + "ACGT"]]


# --------------------------------------------------------------------------- #
# Assignment + tiered overflow
# --------------------------------------------------------------------------- #
def test_assign_identifiers_unique_within_robust_tier():
    codes = ["AAG", "CCA", "GGT", "TTC"]
    ids = finalizeprobes.assign_identifiers(
        [["p1", "p2"], ["q1"]], codes, robust_size=4
    )
    assert ids == [["AAG", "CCA"], ["AAG"]]


def test_assign_identifiers_warns_past_robust_tier():
    codes = ["AAG", "CCA", "GGT", "TTC"]
    with pytest.warns(UserWarning):
        ids = finalizeprobes.assign_identifiers(
            [["p1", "p2", "p3"]], codes, robust_size=2, headers=["GeneX"]
        )
    assert ids == [["AAG", "CCA", "GGT"]]


def test_assign_identifiers_raises_when_over_capacity():
    codes = ["AAG", "CCA"]
    with pytest.raises(ValueError):
        finalizeprobes.assign_identifiers(
            [["p1", "p2", "p3"]], codes, robust_size=2, headers=["GeneX"]
        )


# --------------------------------------------------------------------------- #
# Output column
# --------------------------------------------------------------------------- #
def _write(tmp_path, identifiers):
    dirnames = (str(tmp_path), str(tmp_path / "TempFolder20200101000000"))
    createoutput.writeprobefile(
        ["GeneA"],
        ["GeneA"],
        [["PADLOCK"]],
        [[60.0]],
        [[0]],
        [["TARGET"]],
        dirnames,
        4,
        "out_",
        regions=None,
        identifiers=identifiers,
    )
    out_file = tmp_path / "out_20200101000000.csv"
    return out_file.read_text()


def test_writeprobefile_padlock_id_column_present(tmp_path):
    text = _write(tmp_path, identifiers=[["AAG"]])
    assert "padlock_id" in text.splitlines()[0]
    assert "AAG" in text


def test_writeprobefile_no_padlock_id_column_by_default(tmp_path):
    text = _write(tmp_path, identifiers=None)
    assert "padlock_id" not in text


def _file_path():
    import os

    return os.path.join(os.path.dirname(pi.__file__), "padlock_identifiers_3bp.txt")
