"""E10: can occupancy heterogeneity be separated from switching-rate heterogeneity?

Four ground truths with DIFFERENT mechanisms, two of which share the same clone occupancy distribution:
  OCC   clones differ in occupancy pi_c, shared switching rate k        -> V_pi > 0, V_logk = 0
  RATE  clones share occupancy pi, differ in switching rate k_c         -> V_pi = 0, V_logk > 0
  FROZE clones differ in occupancy, k ~ 0 (stable-state inheritance)    -> V_pi > 0, V_logk = 0, no decay
  BOTH  both vary
OCC and FROZE have the SAME equilibrium occupancy distribution and differ only in switching timescale: the
Milestone 2 plateau cannot tell them apart (Proposition U1), the per-clone estimator should.

Also: a nonstationary organ-induction arm (NONSTAT) where the metastatic compartment applies a common logit shift,
testing whether that shift rescales the plateau (U3).

Writes results/tables/e10_*.csv
"""
import os, sys, time
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.lineage_inference.ctmc import FastTree
from src.lineage_inference.perclone import fit_clone, variance_components, lrt_shared, EXPIT, LOGIT
from src.lineage_inference.propensity import decompose

QUICK = "--quick" in sys.argv
REPS = 3 if QUICK else 8
SEED = 20261128
TAB = os.path.join(ROOT, "results/tables"); os.makedirs(TAB, exist_ok=True)


def sim_clone_tree(pi, k, n_leaves, depth, rng, branch=1.0, miss=0.0):
    """A balanced binary clone tree of the given depth, states evolved down it by a 2-state CTMC."""
    par = [-1]; nodes = [0]
    while len(nodes) < n_leaves and max(len(par) for _ in [0]) < 4 * n_leaves:
        new = []
        for v in nodes:
            for _ in range(2):
                par.append(v); new.append(len(par) - 1)
        nodes = new
        if len(nodes) >= n_leaves:
            break
    par = np.array(par)
    children = {}
    for i, p in enumerate(par):
        children.setdefault(p, []).append(i)
    st = np.full(len(par), -1)
    st[0] = int(rng.random() < pi)
    E = np.exp(-k * branch)
    for i in range(1, len(par)):
        p = st[par[i]]
        up = pi * (1 - E) if p == 0 else pi + (1 - pi) * E
        st[i] = int(rng.random() < up)
    leaves = [i for i in range(len(par)) if i not in children]
    obs = np.full(len(par), -1)
    for v in leaves:
        if rng.random() >= miss:
            obs[v] = st[v]
    return FastTree(par, obs), par, obs


def make(mech, C, rng, n_leaves=64, depth=6):
    trees, truth = [], []
    for c in range(C):
        if mech == "OCC":
            pi, k = float(np.clip(rng.normal(0.5, 0.18), .05, .95)), 0.5
        elif mech == "RATE":
            pi, k = 0.5, float(np.exp(rng.normal(np.log(0.5), 1.0)))
        elif mech == "FROZE":
            pi, k = float(np.clip(rng.normal(0.5, 0.18), .05, .95)), 0.01
        elif mech == "BOTH":
            pi, k = float(np.clip(rng.normal(0.5, 0.18), .05, .95)), float(np.exp(rng.normal(np.log(0.5), 1.0)))
        else:
            raise ValueError(mech)
        t, par, obs = sim_clone_tree(pi, k, n_leaves, depth, rng)
        trees.append(t); truth.append((pi, k))
    return trees, np.array(truth)


def main():
    rng = np.random.default_rng(SEED)
    rows, perclone = [], []
    for mech in ("OCC", "RATE", "FROZE", "BOTH"):
        for C in ((8, 20) if QUICK else (8, 20, 40)):
            reps_c = max(3, REPS // (1 + C // 20))
            for rep in range(reps_c):
                trees, truth = make(mech, C, rng)
                fits = [fit_clone(t) for t in trees]
                vc = variance_components(fits)
                lp = lrt_shared(trees, fits, "pi")
                lk = lrt_shared(trees, fits, "k")
                rows.append(dict(mech=mech, C=C, rep=rep,
                                 true_V_pi=float(np.var(truth[:, 0], ddof=1)),
                                 true_V_logk=float(np.var(np.log(truth[:, 1]), ddof=1)),
                                 lrt_pi_p=lp["p_asym"], lrt_k_p=lk["p_asym"], **vc))
                if rep == 0 and C == 20:
                    for (pt, kt), f in zip(truth, fits):
                        perclone.append(dict(mech=mech, true_pi=pt, true_k=kt, est_pi=f["pi"], est_k=f["k"],
                                             k_at_bound=f["k_at_bound"], n=f["n"]))
            print(mech, C, "done", flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(TAB, "e10_occupancy_vs_rate.csv"), index=False)
    pd.DataFrame(perclone).to_csv(os.path.join(TAB, "e10_perclone_recovery.csv"), index=False)

    # U3: does a common induction shift rescale the plateau?
    ns = []
    for shift in (0.0, 0.5, 1.0, 2.0):
        for rep in range(REPS):
            trees, truth = make("OCC", 20, rng)
            pis = np.clip(EXPIT(LOGIT(np.clip(truth[:, 0], .02, .98)) + shift), .02, .98)
            t2 = []
            for pi, (_, k) in zip(pis, truth):
                t, _, _ = sim_clone_tree(float(pi), float(k), 64, 6, rng)
                t2.append(t)
            fits = [fit_clone(t) for t in t2]
            vc = variance_components(fits)
            ns.append(dict(shift=shift, rep=rep, true_V_pi_shifted=float(np.var(pis, ddof=1)),
                           true_V_pi_base=float(np.var(truth[:, 0], ddof=1)), **vc))
        print("nonstat", shift, flush=True)
    pd.DataFrame(ns).to_csv(os.path.join(TAB, "e10_nonstationary_induction.csv"), index=False)
    print("done")


if __name__ == "__main__":
    t0 = time.time(); main(); print(f"{time.time() - t0:.0f}s")
