"""Per-probe (padlock-specific) identifier sets.

Each padlock for a gene can carry a short identifier so individual probes of the
same gene can be told apart in the readout. An identifier set is an *ordered*
list of equal-length DNA codes with two properties:

1. **Tiered robustness.** The leading ``robust_size`` codes are pairwise Hamming
   distance ``>= min_distance`` (error-robust). The remaining codes are only
   guaranteed distinct (distance >= 1), so a gene can have more padlocks than the
   robust tier with graceful degradation rather than a hard failure.
2. **Informative ordering.** The robust tier is ordered so the first base
   distinguishes as many codes as possible, then the first two bases, and so on.
   This means that for a gene with only a few padlocks you need to sequence only
   the first base or two to tell its probes apart; the trailing bases provide the
   error robustness.

The generation reuses the greedy distance-constrained selection in
:mod:`multi_padlock_design.barcodes.barcodes`.
"""

from itertools import product
from typing import List, Sequence, Tuple

from multi_padlock_design.barcodes import barcodes


def _hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def homopolymer_free(seq: str) -> bool:
    """True if ``seq`` contains no run of 3+ identical consecutive bases."""
    return not any(b * 3 in seq for b in "ACGT")


def candidate_pool(length: int, exclude_homopolymers: bool = True) -> List[str]:
    """All ``length``-mers over ACGT, optionally excluding homopolymer runs."""
    pool = ["".join(p) for p in product("ACGT", repeat=length)]
    if exclude_homopolymers:
        pool = [s for s in pool if homopolymer_free(s)]
    return pool


def informative_order(codes: Sequence[str]) -> List[str]:
    """Reorder ``codes`` so short prefixes are maximally discriminative.

    Greedily picks, at each step, the code that collides least with the
    already-chosen codes at the shortest prefix first (fewest sharing the same
    first base, then the same two bases, ...). This cycles the first base as fast
    as possible so the first block of codes is separable on base 1 alone, the
    next prefix length separates the next tier, etc. Ties break alphabetically
    for determinism.
    """
    codes = list(codes)
    if not codes:
        return []
    length = len(codes[0])
    chosen: List[str] = []
    remaining = list(codes)
    while remaining:

        def collision_key(c: str):
            prefix_collisions = tuple(
                sum(1 for o in chosen if o[:p] == c[:p]) for p in range(1, length + 1)
            )
            return prefix_collisions + (c,)  # alphabetical tie-break

        best = min(remaining, key=collision_key)
        chosen.append(best)
        remaining.remove(best)
    return chosen


def extend_distinct(robust: Sequence[str], pool: Sequence[str]) -> List[str]:
    """Append the remaining ``pool`` codes, most-spread first.

    Each appended code is the one maximizing its minimum Hamming distance to the
    codes chosen so far. Only distinctness (distance >= 1) is guaranteed for the
    extension; it simply keeps the extra codes as separated as possible.
    """
    chosen = list(robust)
    chosen_set = set(chosen)
    remaining = [c for c in pool if c not in chosen_set]
    while remaining:
        best = max(remaining, key=lambda c: min(_hamming(c, x) for x in chosen))
        chosen.append(best)
        remaining.remove(best)
    return chosen


def build_identifier_list(
    length: int = 3,
    min_distance: int = 2,
    exclude_homopolymers: bool = True,
    restarts: int = 500,
    seed: int = 7,
) -> Tuple[List[str], int]:
    """Build an ordered, tiered identifier set.

    Returns ``(codes, robust_size)`` where ``codes[:robust_size]`` is the
    informatively-ordered robust tier (pairwise distance >= ``min_distance``) and
    the remainder is the graceful distinct extension.
    """
    pool = candidate_pool(length, exclude_homopolymers=exclude_homopolymers)
    robust = barcodes.multi_restart_greedy(
        n=length,
        alphabet="ACGT",
        d=min_distance,
        restarts=restarts,
        seed=seed,
        candidates=list(pool),
    )
    robust = informative_order(robust)
    codes = extend_distinct(robust, pool)
    return codes, len(robust)


def load_identifiers(path: str) -> List[str]:
    """Read one code per line (matching the existing barcode-file format)."""
    with open(path) as f:
        return [line.strip().upper() for line in f if line.strip()]


def robust_prefix_size(codes: Sequence[str], min_distance: int) -> int:
    """Length of the leading prefix of ``codes`` pairwise >= ``min_distance``."""
    for k, ck in enumerate(codes):
        if any(_hamming(ck, codes[j]) < min_distance for j in range(k)):
            return k
    return len(codes)


def validate_identifier_set(codes: Sequence[str], min_distance: int) -> int:
    """Validate an identifier set and return its robust-prefix size.

    Raises ``ValueError`` if the set is empty, codes have unequal length, or
    contains duplicates.
    """
    if not codes:
        raise ValueError("Identifier set is empty.")
    lengths = {len(c) for c in codes}
    if len(lengths) != 1:
        raise ValueError(f"Identifier codes have unequal lengths: {sorted(lengths)}")
    if len(set(codes)) != len(codes):
        raise ValueError("Identifier set contains duplicate codes.")
    return robust_prefix_size(codes, min_distance)


def write_identifiers(path: str, codes: Sequence[str]) -> None:
    """Write codes one per line (trailing newline)."""
    with open(path, "w") as f:
        f.write("\n".join(codes) + "\n")


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(
        description="Generate a tiered padlock identifier set."
    )
    ap.add_argument("--length", type=int, default=3)
    ap.add_argument("--min-distance", type=int, default=2)
    ap.add_argument("--restarts", type=int, default=500)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--allow-homopolymers", action="store_true")
    ap.add_argument(
        "--out", required=True, help="Output .txt path (one code per line)."
    )
    args = ap.parse_args()

    codes, robust = build_identifier_list(
        length=args.length,
        min_distance=args.min_distance,
        exclude_homopolymers=not args.allow_homopolymers,
        restarts=args.restarts,
        seed=args.seed,
    )
    write_identifiers(args.out, codes)
    print(
        f"Wrote {len(codes)} identifiers to {args.out} "
        f"(robust tier of {robust} at min distance {args.min_distance})."
    )


if __name__ == "__main__":
    main()
