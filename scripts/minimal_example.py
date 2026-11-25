"""Minimal working example: simulate clones that differ only in occupancy, then recover it.

Twelve clones share a switching rate and differ in how much of the time they hold the
state. The estimator should recover the occupancy spread and should NOT report a
meaningful spread of rates.

Run from the repository root:

    python -m scripts.minimal_example
    python scripts/minimal_example.py --clones 12 --cells 48 --seed 0
"""
import argparse
import numpy as np

from src.branching.recorder import simulate_genealogy, induced_subtree
from src.lineage_inference.ctmc import VarTree, fit_var

K_SHARED = 0.3          # every clone switches at the same rate, so true V_log k = 0


def evolve(parent, rng, pi, k):
    """Two-state CTMC down the genealogy, one unit of time per generation."""
    E = np.exp(-k)
    st = np.empty(len(parent), int)
    st[0] = int(rng.random() < pi)
    for i in range(1, len(parent)):
        p = st[parent[i]]
        up = pi * (1 - E) if p == 0 else pi + (1 - pi) * E
        st[i] = int(rng.random() < up)
    return st


def one_clone(rng, pi, n_gen, n_cells):
    par, gen, live = simulate_genealogy(n_gen, rng)
    n_cells = min(n_cells, len(live))
    take = rng.choice(live, n_cells, replace=False)
    sub_par, sub_elen, leafmap, _ = induced_subtree(par, take)
    st = evolve(par, rng, pi, K_SHARED)
    leaf = np.full(len(sub_par), -1, int)
    for new_i, orig in leafmap.items():
        leaf[new_i] = st[orig]
    return VarTree(sub_par, leaf, elen=sub_elen)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clones", type=int, default=12)
    ap.add_argument("--cells", type=int, default=48, help="cells sampled per clone")
    ap.add_argument("--generations", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    true_pi = rng.uniform(0.15, 0.85, a.clones)

    print(f"simulating {a.clones} clones x {a.cells} cells, {a.generations} generations")
    print(f"every clone shares k = {K_SHARED}, so the true across-clone rate variance is 0\n")
    print(f"{'clone':>6} {'true pi':>9} {'est pi':>8} {'est k':>8}  {'at bound':>9}")

    est_pi, est_k, flagged = [], [], 0
    for c in range(a.clones):
        fit = fit_var(one_clone(rng, true_pi[c], a.generations, a.cells))
        est_pi.append(fit["pi"])
        est_k.append(fit["k"])
        flagged += int(fit["k_at_bound"])
        print(f"{c:>6} {true_pi[c]:>9.3f} {fit['pi']:>8.3f} {fit['k']:>8.3f}"
              f"  {'yes' if fit['k_at_bound'] else '':>9}")

    est_pi = np.asarray(est_pi)
    est_k = np.asarray(est_k)
    print("\n" + "-" * 58)
    print(f"{'V_pi   (occupancy variance)':<34} true {np.var(true_pi):.4f}"
          f"   est {np.var(est_pi):.4f}")
    print(f"{'V_log k (rate variance)':<34} true {0.0:.4f}"
          f"   est {np.var(np.log(est_k)):.4f}")
    print(f"{'correlation of est pi with truth':<34} {np.corrcoef(est_pi, true_pi)[0,1]:.3f}")
    print(f"{'clones with k at the bound':<34} {flagged}/{a.clones}")
    print("-" * 58)
    print("\nOccupancy is recovered. The estimated rate variance is NOT zero even though the")
    print("truth is: a spread of per-clone rates is not evidence that clones differ in rate.")
    print("See docs/methodology.md and benchmarks/reports/benchmark_results.md.")


if __name__ == "__main__":
    main()
