"""E11: optimal MRCA-depth sampling (Proposition U5), analytic Fisher information plus an equal-budget simulation.

Analytic part: the Fisher information a single within-clone pair at lineage distance L carries about the switching
rate k, and about the occupancy pi. Simulation part: equal total-cell budgets spent in different ways.
"""
import os, sys, time
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.lineage_inference.perclone import fit_clone, variance_components

TAB = os.path.join(ROOT, "results/tables"); os.makedirs(TAB, exist_ok=True)
QUICK = "--quick" in sys.argv


def pair_probs(pi, rho):
    q = 1 - pi
    p11 = pi * pi + pi * q * rho
    p00 = q * q + pi * q * rho
    p10 = p01 = pi * q * (1 - rho)
    return np.array([p00, p01, p10, p11])


def info_pair(pi, k, L):
    """Fisher information of one binary pair about (log k, pi). Returns (I_logk, I_pi)."""
    q = 1 - pi
    rho = np.exp(-2 * k * L)
    p = pair_probs(pi, rho)
    p = np.clip(p, 1e-12, None)
    dp_drho = np.array([pi * q, -pi * q, -pi * q, pi * q])
    I_rho = float((dp_drho ** 2 / p).sum())
    I_logk = I_rho * (2 * k * L * rho) ** 2                      # drho/dlogk = -2kL rho
    dp_dpi = np.array([-2 * q + (1 - 2 * pi) * rho, (1 - 2 * pi) * (1 - rho),
                       (1 - 2 * pi) * (1 - rho), 2 * pi + (1 - 2 * pi) * rho])
    I_pi = float((dp_dpi ** 2 / p).sum())
    return I_logk, I_pi


def main():
    rows = []
    for k in (0.1, 0.5, 2.0):
        for pi in (0.3, 0.5):
            for L in np.exp(np.linspace(np.log(0.02), np.log(60), 120)):
                a, b = info_pair(pi, k, L)
                rows.append(dict(k=k, pi=pi, L=L, kL=k * L, I_logk=a, I_pi=b))
    d = pd.DataFrame(rows); d.to_csv(os.path.join(TAB, "e11_fisher_information.csv"), index=False)
    opt = d.loc[d.groupby(["k", "pi"]).I_logk.idxmax(), ["k", "pi", "L", "kL", "I_logk"]]
    opt["L_star_theory"] = 1 / (2 * opt.k)
    opt.to_csv(os.path.join(TAB, "e11_optimal_depth.csv"), index=False)
    print(opt.round(4).to_string(), flush=True)

    # equal-budget designs: the same total cells spent on clones x cells x depth spread
    rng = np.random.default_rng(20261129)
    from experiments.lineage_identifiability.run_e10 import sim_clone_tree
    reps = 3 if QUICK else 10
    budget = 1280
    designs = [("many_clones_shallow", 40, 32, 5), ("balanced", 20, 64, 6),
               ("few_clones_deep", 10, 128, 7), ("very_few_very_deep", 5, 256, 8)]
    out = []
    for name, C, ncell, depth in designs:
        for rep in range(reps):
            trees, truth = [], []
            for c in range(C):
                pi = float(np.clip(rng.normal(0.5, 0.18), .05, .95))
                kk = float(np.exp(rng.normal(np.log(0.5), 1.0)))
                t, _, _ = sim_clone_tree(pi, kk, ncell, depth, rng)
                trees.append(t); truth.append((pi, kk))
            truth = np.array(truth)
            fits = [fit_clone(t) for t in trees]
            vc = variance_components(fits)
            out.append(dict(design=name, C=C, cells_per_clone=ncell, depth=depth, total=C * ncell, rep=rep,
                            true_V_pi=float(np.var(truth[:, 0], ddof=1)),
                            true_V_logk=float(np.var(np.log(truth[:, 1]), ddof=1)), **vc))
        print(name, "done", flush=True)
    pd.DataFrame(out).to_csv(os.path.join(TAB, "e11_equal_budget.csv"), index=False)
    print("done")


if __name__ == "__main__":
    t0 = time.time(); main(); print(f"{time.time() - t0:.0f}s")
