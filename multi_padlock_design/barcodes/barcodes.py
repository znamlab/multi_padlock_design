#!/usr/bin/env python3
import argparse
import multiprocessing as mp
import os
import random
import time
from itertools import product

import numpy as np

BASE2INT = {'A':0, 'C':1, 'G':2, 'T':3}
INT2BASE = np.array(['A','C','G','T'])

# ---- Globals to avoid re-sending big data to workers ----
_GLOBAL_ENCODED = None

## Helpers
def encode_str(s: str) -> np.ndarray:
    return np.fromiter((BASE2INT[ch] for ch in s), dtype=np.uint8)

def encode_list(strs) -> np.ndarray:
    if len(strs) == 0:
        return np.empty((0,0), dtype=np.uint8)
    encoded = []
    for s in strs:
        if isinstance(s, str):
            encoded.append(encode_str(s))
        else:
            encoded.append(np.asarray(s, dtype=np.uint8))
    return np.vstack(encoded)

def decode_arr(arr: np.ndarray) -> str:
    arr = np.asarray(arr, dtype=np.uint8)   # make sure it’s an array
    return ''.join(INT2BASE[arr.tolist()])  # convert safely

def hamming_vectorized(X, y):
    return (X != y).sum(axis=1)

def hamming_dists_to_one(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    # X: (m, n), y: (n,)
    if X.size == 0:
        return np.empty((0,), dtype=np.int32)
    return (X != y).sum(axis=1)

def hamming_distance(a: np.ndarray, b: np.ndarray) -> int:
    return np.count_nonzero(a != b)

## Library constrains

def make_constrained_library(repeat):
    #filtering
    #Ensure there are no runs of more than `max_homopolymer_length` consecutive bases
    #and at least three different bases

    #look at all 10_mers
    candidates = np.array([''.join(p) for p in product("ACGT", repeat=repeat)])
    print(f'Length before constraints {len(candidates)}')
    #homopolymer constraint
    def noHomo(barcode):
        if ('AAA' in barcode) | ('TTT' in barcode) | ('CCC' in barcode) | ("GGG" in barcode):
            noHomo=False
        else:
            noHomo=True
        return noHomo

    noHomo_mask = [noHomo(barcode) for barcode in candidates]
    candidates = candidates[noHomo_mask]

    print(f'length after homopolymers: {len(candidates)}')

    diversity_mask = [(len(set(barcode))>=3) for barcode in candidates]
    candidates = candidates[diversity_mask]

    print(f'length after diversity: {len(candidates)}')

    return candidates

## Build first greedy panel

def greedy_once(candidates, d):
    chosen = []
    X = _GLOBAL_ENCODED
    for i in idxs:
        cand = X[i]
        if all(hamming_distance(cand, c) >= d for c in chosen):
            chosen.append(cand)
    return chosen

def multi_restart_greedy(n=7, alphabet="ACGT", d=3, restarts=20, seed=None, candidates=None):
    #really need to make this parallelizable, capitalise on the cluster ffs
    print(f'multi-restart greedy: {restarts} iterations')
    if candidates is None:
        print("no candidates were supplied, building complete library")
        candidates = [encode_str(''.join(p)) for p in product(alphabet, repeat=n)]
    else:
        encoded = encode_list(candidates)

    best = []
    rng = random.Random(seed)
    start = time.time()
    for _ in range(restarts):
        print(f"Iteration{_}, {start-time.time():.2f} s/it")
        rng.shuffle(candidates)
        chosen = greedy_once(candidates, d)
        if len(chosen) > len(best):
            best = chosen
        start = time.time()
    # decode back if you want strings
    return [decode_arr(c) for c in best]


def _init_candidates(encoded):
    # store a reference in each worker (cheap under fork; acceptable under spawn)
    global _GLOBAL_ENCODED
    _GLOBAL_ENCODED = encoded

def _worker_greedy(args):
    print(f"Worker starting seed {args[-1]}")
    # args: (d, seed)
    d, seed = args
    rng = random.Random(seed)
    # shuffle indices instead of the objects
    idxs = list(range(len(_GLOBAL_ENCODED)))
    rng.shuffle(idxs)
    # greedy over the shuffled view without copying candidates
    chosen = []
    for i in idxs:
        cand = _GLOBAL_ENCODED[i]
        if all(hamming_distance(cand, c) >= d for c in chosen):
            chosen.append(cand)
    return len(chosen), chosen  # length first to reduce parent work

def parallel_greedy(n=7, alphabet="ACGT", d=3, restarts=20, candidates=None, base_seed=123):
    """
    Run greedy in parallel across multiple random shuffles of the candidate list.
    Returns: (best_decoded_list, distribution_of_sizes)
    """
    # Prepare encoded candidates once in parent
    if candidates is None:
        # Full 4^n is gigantic when n=13; I assume your constrained builder is used.
        candidates = [''.join(p) for p in product(alphabet, repeat=n)]
    encoded = [encode_str(c) for c in candidates]

    # Pool size aligned to SLURM allocation
    slurm_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) or os.cpu_count() or 1
    pool_size = min(restarts, slurm_cpus)

    seeds = [base_seed + i for i in range(restarts)]
    worker_args = [(d, s) for s in seeds]

    # Prefer 'fork' on Linux to avoid copying big data; spawn is safer but slower
    # mp.set_start_method('fork', force=True)  # optional, Linux only

    t0 = time.time()
    print(f"multi-restart greedy: {restarts} iterations (pool={pool_size})")

    with mp.get_context("fork").Pool(
        processes=pool_size,
        initializer=_init_candidates,
        initargs=(encoded,),
        maxtasksperchild=1  # keeps memory tidy on long runs
    ) as pool:
        results = pool.map(_worker_greedy, worker_args)

    # pick best by length
    lengths = [r[0] for r in results]
    best_len, best_choice = max(results, key=lambda x: x[0])
    distribution = lengths

    print(f"Finished in {time.time()-t0:.2f}s. Best={best_len}, mean={sum(lengths)/len(lengths):.1f}")

    # decode only the winning set
    best_decoded = [decode_arr(c) for c in best_choice]
    return best_decoded, distribution


## Hill climbing

def local_improvement_fast(chosen, candidates, d, max_passes=3):
    print(f'Local improvement(hill-climbing). Max: {max_passes} passes')
    # Encode chosen and candidates
    chosen_enc = encode_list(chosen)
    cand_enc   = encode_list(candidates)
    chosen_set = set(chosen)  # quick membership check
    start = time.time()

    for _ in range(max_passes):
        print(f"Pass {_}, at {time.time()-start} s/pass")
        improved = False
        new_chosen = []
        for x_str, x_enc in zip(candidates, cand_enc):
            if x_str in chosen_set:
                continue
            # vectorized distance check
            dists = (chosen_enc != x_enc).sum(axis=1)
            conflicts = np.where(dists < d)[0]

            if len(conflicts) == 0:
                chosen.append(x_str)
                chosen_enc = np.vstack([chosen_enc, x_enc])
                chosen_set.add(x_str)
                improved = True
            elif len(conflicts) == 1:
                # attempt swap
                idx = conflicts[0]
                tmp_enc = np.delete(chosen_enc, idx, axis=0)
                if ((tmp_enc != x_enc).sum(axis=1) >= d).all():
                    # swap is valid
                    removed = chosen[idx]
                    chosen.pop(idx)
                    chosen_enc = tmp_enc
                    chosen.append(x_str)
                    chosen_enc = np.vstack([chosen_enc, x_enc])
                    chosen_set.discard(removed)
                    chosen_set.add(x_str)
                    improved = True
        if not improved:
            break

        start = time.time()

    return chosen

def worker_local_improvement(args):
    chosen, candidates, d, max_passes, seed = args
    rng = random.Random(seed)
    # shuffle candidates independently for each worker
    shuffled_candidates = list(candidates)
    rng.shuffle(shuffled_candidates)
    improved = local_improvement_fast(list(chosen), shuffled_candidates, d, max_passes=max_passes)
    return improved

def parallel_local_improvement(chosen, candidates, d, max_passes=3, processes=4, base_seed=123):
    """
    Run local_improvement_fast in parallel with multiple random shuffles of candidates.
    Returns the best improved set (largest size).
    """
    seeds = [base_seed + i for i in range(processes)]
    worker_args = [(chosen, candidates, d, max_passes, seed) for seed in seeds]

    with mp.get_context("spawn").Pool(processes=processes) as pool:
        results = pool.map(worker_local_improvement, worker_args)

    # pick the largest improved set
    best = max(results, key=len)
    distribution = [len(result) for result in results]

    return best, distribution

# parallelization

# ------------------- Subset selection (prefix 7, d>=3) -------------------
def maximize_prefix_subset(all_code, d_prefix=3, restarts=16, seed=123):
    """Return the largest subset (greedy w/ restarts) whose first 7 bases
       are pairwise Hamming distance >= d_prefix."""
    rng = random.Random(seed)
    pref7 = np.vstack([encode_str(s[:7]) for s in all_code])  # (m,7)
    idxs = list(range(len(all_code)))
    best = []
    for _ in range(restarts):
        rng.shuffle(idxs)
        chosen_idx = []
        chosen_enc = np.empty((0,7), dtype=np.uint8)
        for i in idxs:
            y = pref7[i]
            if chosen_enc.size == 0:
                chosen_idx.append(i)
                chosen_enc = y.reshape(1,7)
            else:
                # accept if min distance to chosen >= d_prefix
                if ( (chosen_enc != y).sum(axis=1) >= d_prefix ).all():
                    chosen_idx.append(i)
                    chosen_enc = np.vstack([chosen_enc, y])
        if len(chosen_idx) > len(best):
            best = chosen_idx
    return [all_code[i] for i in best]

# ------------------- Simple hill-climber on global set -------------------
def local_improve_simple_fast(all_code,
                              d_global=4, d_prefix=3, alpha=3.0,
                              iters=3000, seed=123,
                              subset_restarts=8, constraints=False):
    """
    Modifies all_code (10-mers) by random insertions or 1-for-1 swaps.
    Accepts a change only if it increases score = alpha*|subset| + |all|.
    If constraints=True, candidates are drawn from constrained library.
    Otherwise, propose fully random 10-mers.
    """
    rng = random.Random(seed)

    # prepare candidate pool if constrained
    cand_pool = None
    if constraints:
        cand_pool = make_constrained_library(10)
        print(f"Constrained pool size: {len(cand_pool)}")

    # pre-encode and score
    all_enc = encode_list(all_code)
    all_set = set(all_code)

    subset = maximize_prefix_subset(all_code, d_prefix=d_prefix,
                                    restarts=subset_restarts, seed=seed)
    cur_score = alpha * len(subset) + len(all_code)

    seen = set()

    for _ in range(iters):
        # --- candidate proposal ---
        if constraints:
            s = rng.choice(cand_pool)
            y = encode_str(s)
        else:
            y = np.fromiter((rng.randrange(4) for _ in range(10)), dtype=np.uint8)
            s = ''.join('ACGT'[b] for b in y)

        if s in all_set or s in seen:
            continue
        seen.add(s)

        dists = (all_enc != y).sum(axis=1)
        conflicts = np.where(dists < d_global)[0]

        # Try insertion (no conflicts)
        if conflicts.size == 0:
            trial_all = all_code + [s]
            trial_subset = maximize_prefix_subset(
                trial_all, d_prefix=d_prefix,
                restarts=subset_restarts, seed=rng.randrange(10**9)
            )
            new_score = alpha * len(trial_subset) + len(trial_all)
            if new_score > cur_score:
                # accept
                all_code = trial_all
                all_enc = np.vstack([all_enc, y])
                all_set.add(s)
                subset = trial_subset
                cur_score = new_score
            continue

        # Try simple 1-for-1 swap (exactly one conflict)
        if conflicts.size == 1:
            victim_idx = conflicts[0]
            victim = all_code[victim_idx]
            trial_all = all_code[:victim_idx] + [s] + all_code[victim_idx+1:]
            trial_subset = maximize_prefix_subset(
                trial_all, d_prefix=d_prefix,
                restarts=subset_restarts, seed=rng.randrange(10**9)
            )
            new_score = alpha * len(trial_subset) + len(trial_all)
            if new_score > cur_score:
                # accept swap
                all_code = trial_all
                all_enc[victim_idx, :] = y
                all_set.remove(victim)
                all_set.add(s)
                subset = trial_subset
                cur_score = new_score

    return all_code, subset, cur_score
# ------------------- Worker + Orchestration -------------------
def worker_run(args_tuple):
    """One independent hill-climb with its own seed; returns (score, all, subset, seed)."""
    (seed, all_code_init, d_global, d_prefix, alpha, iters, subset_restarts) = args_tuple

    # Avoid OpenMP/MKL oversubscription inside workers (one thread per proc)
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    # Copy the starting set so each worker starts identically
    all_code = list(all_code_init)
    improved_all, improved_subset, final_score = local_improve_simple_fast(
        all_code,
        d_global=d_global,
        d_prefix=d_prefix,
        alpha=alpha,
        iters=iters,
        seed=seed,
        subset_restarts=subset_restarts
    )
    return (final_score, improved_all, improved_subset, seed)

def main():
    ap = argparse.ArgumentParser(description="Parallel nested barcode optimizer (8 processes, 8 seeds).")
    ap.add_argument("--input", required=True, help="Path to initial 10-mer list (CSV or txt). "
                                                   "If CSV (id,seq), provide --col to pick the sequence column.")
    ap.add_argument("--col", default=None, help="Column name or index for sequence if input is CSV.")
    ap.add_argument("--delimiter", default=",", help="Delimiter for CSV (default ',').")
    ap.add_argument("--outdir", default="nested_barcodes_parallel_out", help="Output directory.")
    ap.add_argument("--processes", type=int, default=8, help="Number of parallel processes (default 8).")
    ap.add_argument("--base-seed", type=int, default=42, help="Base seed (each worker uses base+idx).")
    ap.add_argument("--iters", type=int, default=5000, help="Iterations per worker.")
    ap.add_argument("--subset-restarts", type=int, default=8, help="Greedy restarts when computing subset.")
    ap.add_argument("--alpha", type=float, default=3.0, help="Score weight: alpha*|subset| + |all|.")
    ap.add_argument("--d-global", type=int, default=4, help="Global Hamming distance constraint (10-mers).")
    ap.add_argument("--d-prefix", type=int, default=3, help="Prefix (first 7) min distance for subset.")
    ap.add_argument("--is-csv", action="store_true", help="Treat input as CSV; otherwise read one sequence per line "
                                                          "(or id,seq per line and take the second field).")
    args = ap.parse_args()

    # Load initial all_code
    all_code = []
    if args.is_csv:
        import csv
        col = args.col
        with open(args.input, newline="") as f:
            reader = csv.reader(f, delimiter=args.delimiter)
            for row in reader:
                if not row:
                    continue
                if isinstance(col, str):
                    # if header row exists, you should adapt to DictReader; for now assume index
                    raise ValueError("--col as name needs DictReader; either provide index or use a TSV with known index.")
                idx = int(col) if col is not None else 1  # default to second column
                seq = row[idx].strip()
                if seq:
                    all_code.append(seq)
    else:
        # text: allow "id,seq" or bare "seq"
        with open(args.input) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if args.delimiter in line:
                    parts = line.split(args.delimiter)
                    seq = parts[-1].strip()
                else:
                    seq = line
                if seq:
                    all_code.append(seq)

    if not all_code:
        raise ValueError("No sequences loaded from --input")

    # Prepare output dir
    os.makedirs(args.outdir, exist_ok=True)

    # Build worker argument tuples
    seeds = [args.base_seed + i for i in range(args.processes)]
    worker_args = [
        (seed, all_code, args.d_global, args.d_prefix, args.alpha, args.iters, args.subset_restarts)
        for seed in seeds
    ]

    # Run in parallel
    with mp.get_context("spawn").Pool(processes=args.processes) as pool:
        results = pool.map(worker_run, worker_args)

    # Pick the best
    results.sort(key=lambda x: x[0], reverse=True)  # by score desc
    best_score, best_all, best_subset, best_seed = results[0]

    # Save
    all_path = os.path.join(args.outdir, f"{len(best_all)}_d{args.d_global}_10bp_barcodes_seed{best_seed}.txt")
    subset_path = os.path.join(args.outdir, f"{len(best_subset)}_subset_prefix7_d{args.d_prefix}_seed{best_seed}.txt")
    summary_path = os.path.join(args.outdir, "summary.txt")

    with open(all_path, "w") as f:
        for i, code in enumerate(best_all):
            f.write(f"{i},{code}\n")

    with open(subset_path, "w") as f:
        for i, code in enumerate(best_subset):
            f.write(f"{i},{code}\n")

    with open(summary_path, "w") as f:
        f.write(f"Processes: {args.processes}\n")
        f.write(f"Seeds: {seeds}\n")
        f.write(f"Best seed: {best_seed}\n")
        f.write(f"Global |all|: {len(best_all)}\n")
        f.write(f"Subset |subset|: {len(best_subset)}\n")
        f.write(f"Score (alpha*|subset| + |all|): {best_score}\n")
        f.write(f"d_global: {args.d_global}, d_prefix: {args.d_prefix}, alpha: {args.alpha}\n")
        f.write(f"iters per worker: {args.iters}, subset_restarts: {args.subset_restarts}\n")

    print("Saved:")
    print("  Global set:", all_path)
    print("  Subset:    ", subset_path)
    print("  Summary:   ", summary_path)