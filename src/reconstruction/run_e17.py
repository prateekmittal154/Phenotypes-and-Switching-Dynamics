"""E17: where is switching information lost? Decomposition across the T0-T1-T2-T3 hierarchy.

For each replicate a full genealogy (T0) is simulated with states and barcodes, cells are sampled (T1), the
recorder output is produced (T2) and a tree is reconstructed (T3). The switching rate is then estimated four
ways on the SAME observed leaf states:

  L1  true sampled topology with TRUE branch lengths        (ancestry perfectly known)
  L2  true sampled topology with UNIT edges                 (adds branch-length miscalibration only)
  L3  reconstructed topology with unit edges                (adds topology error: the analyst's situation)
  L3e ensemble over reconstructions (greedy / nj / upgma), averaging the per-tree likelihood

Comparing L1 to L2 isolates branch-length miscalibration; L2 to L3 isolates topology error; and the spread of
k-hat gives an operational information measure, since Fisher information is the inverse of the estimator's
variance for an efficient estimator.

Usage: python -I run_e17.py [--quick]
"""
import os, sys, time
from multiprocessing import Pool
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.branching.recorder import (simulate_genealogy, evolve_states, induced_subtree, record_barcodes,
                                    reconstruct, mrca_depths, colless, saturation, MISSING)
from src.lineage_inference.ctmc import VarTree, fit_var

TAB = os.path.join(ROOT, "results/tables"); os.makedirs(TAB, exist_ok=True)
QUICK = "--quick" in sys.argv
REPS = 6 if QUICK else 40
N_SAMPLE = 96
N_GEN = 10

# Two recorder settings. Provenance and uncertainty are documented in
# reports/milestone5_recorder_simulation.md; these are plausible stand-ins, not fitted to the real data.
RECORDERS = {
    "gestalt_like": dict(n_sites=9, mu=0.30, n_states=12, dropout=0.01, missing=0.12),
    "cassiopeia_like": dict(n_sites=25, mu=0.12, n_states=40, dropout=0.005, missing=0.20),
}
GENEALOGIES = {
    "balanced": dict(p_death=0.02, bottleneck_at=None, bottleneck_size=None),
    "bottlenecked": dict(p_death=0.02, bottleneck_at=4, bottleneck_size=4),
}


def _clades(par, n_leaves):
    """Set of non-trivial leaf-sets defined by internal nodes, for Robinson-Foulds style comparison."""
    par = np.asarray(par)
    ch = {}
    for v, p in enumerate(par):
        if p >= 0:
            ch.setdefault(p, []).append(v)
    out = set()

    def rec(v):
        if v not in ch:
            return frozenset([v]) if v < n_leaves else frozenset()
        s = frozenset().union(*[rec(u) for u in ch[v]])
        if 1 < len(s) < n_leaves:
            out.add(s)
        return s

    root = int(np.where(par < 0)[0][0])
    rec(root)
    return out


def one(args):
    rec_name, gen_name, k_true, rep = args
    rng = np.random.default_rng(abs(hash((rec_name, gen_name, round(k_true, 4), rep))) % (2**32))
    rp = RECORDERS[rec_name]; gp = GENEALOGIES[gen_name]
    pi_true = 0.4
    par, gen, live = simulate_genealogy(N_GEN, rng, **gp)
    if len(live) < N_SAMPLE:
        return None
    st = evolve_states(par, rng, pi_true, k_true)
    samp = rng.choice(live, N_SAMPLE, replace=False)
    # ---- T1
    p1, l1, leafmap, order = induced_subtree(par, samp)
    leaf_state = np.full(len(p1), -1)
    for newi, orig in leafmap.items():
        leaf_state[newi] = st[orig]
    obs = np.array([st[s] for s in samp])
    # ---- T2
    A = record_barcodes(par, gen, rng, **rp)
    Al = A[samp]
    sat = saturation(Al)
    n_informative = int(sum(1 for c in range(Al.shape[1])
                            for s in np.unique(Al[:, c][Al[:, c] > 0])
                            if 0 < (Al[:, c] == s).sum() < len(Al)))
    # ---- fits
    out = dict(recorder=rec_name, genealogy=gen_name, k_true=k_true, rep=rep,
               saturation=sat, missing=float((Al == MISSING).mean()), n_informative_chars=n_informative,
               n_leaves=N_SAMPLE)
    t1 = VarTree(p1, leaf_state, elen=l1)
    f1 = fit_var(t1); out.update(k_L1=f1["k"], pi_L1=f1["pi"], bound_L1=f1["k_at_bound"])
    t2 = VarTree(p1, leaf_state, elen=np.ones(len(p1)))
    f2 = fit_var(t2); out.update(k_L2=f2["k"], pi_L2=f2["pi"], bound_L2=f2["k_at_bound"])
    d1, paths1 = mrca_depths(p1, len(leafmap), rng=rng)
    h1 = max(len(q) for q in paths1) - 1
    out.update(mrca_mean_T1=float(d1.mean()), mrca_norm_T1=float(d1.mean() / max(h1, 1)),
               height_T1=h1, colless_T1=colless(p1, len(leafmap)),
               close_frac_T1=float((d1 >= 0.7 * h1).mean()))
    clades1 = _clades(p1, len(leafmap))
    lls = {}
    for method in ("greedy", "nj", "upgma"):
        p3 = reconstruct(Al, method, rng)
        ls3 = np.full(len(p3), -1); ls3[:N_SAMPLE] = obs
        t3 = VarTree(p3, ls3)
        f3 = fit_var(t3)
        out[f"k_L3_{method}"] = f3["k"]; out[f"bound_L3_{method}"] = f3["k_at_bound"]
        d3, paths3 = mrca_depths(p3, N_SAMPLE, rng=rng)
        h3 = max(len(q) for q in paths3) - 1
        out[f"mrca_mean_{method}"] = float(d3.mean())
        out[f"mrca_norm_{method}"] = float(d3.mean() / max(h3, 1))
        out[f"height_{method}"] = h3
        out[f"colless_{method}"] = colless(p3, N_SAMPLE)
        out[f"close_frac_{method}"] = float((d3 >= 0.7 * h3).mean())
        c3 = _clades(p3, N_SAMPLE)
        inter = len(clades1 & c3)
        out[f"rf_recall_{method}"] = inter / max(len(clades1), 1)
        out[f"rf_precision_{method}"] = inter / max(len(c3), 1)
        lls[method] = (t3, f3["pi"])
    # ---- ensemble: average the likelihood over the three reconstructions
    def ens_nll(pi, k):
        v = [t.loglik(pi, pi, k) for t, _ in lls.values()]
        m = max(v)
        return -(m + np.log(np.mean(np.exp(np.array(v) - m))))
    best = (np.inf, np.nan, np.nan)
    for p in np.linspace(0.08, 0.92, 15):
        for k in np.exp(np.linspace(np.log(5e-3), np.log(15), 20)):
            v = ens_nll(p, k)
            if v < best[0]:
                best = (v, p, k)
    out["k_L3_ens"] = best[2]; out["pi_L3_ens"] = best[1]
    return out


def main():
    jobs = [(r, g, k, rep) for r in RECORDERS for g in GENEALOGIES
            for k in (0.3, 1.0, 2.0) for rep in range(REPS)]
    with Pool(8) as p:
        rows = [x for x in p.map(one, jobs) if x]
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(TAB, "e17_hierarchy.csv"), index=False)

    def summ(g, kt_val):
        o = {}
        for lab in ("L1", "L2", "L3_greedy", "L3_nj", "L3_upgma", "L3_ens"):
            col = f"k_{lab}"
            if col not in g:
                continue
            kk = g[col].values; kt = kt_val
            o[f"{lab}_median"] = float(np.median(kk))
            o[f"{lab}_logbias"] = float(np.median(np.log(kk / kt)))
            o[f"{lab}_logsd"] = float(np.std(np.log(kk)))
        o["n"] = len(g)
        return pd.Series(o)

    recs = []
    for (r_, g_, kt_), grp in d.groupby(["recorder", "genealogy", "k_true"]):
        row = summ(grp, kt_)
        row["recorder"], row["genealogy"], row["k_true"] = r_, g_, kt_
        recs.append(row)
    s = pd.DataFrame(recs)
    s.to_csv(os.path.join(TAB, "e17_hierarchy_summary.csv"), index=False)
    pd.set_option("display.width", 250)
    print(s.round(3).to_string(index=False))
    print("\nrelative information (Var of log k-hat at L1) / (Var at level), >1 means information lost:")
    for lab in ("L2", "L3_greedy", "L3_nj", "L3_upgma", "L3_ens"):
        r = (s[f"L1_logsd"] ** 2 / s[f"{lab}_logsd"] ** 2).median()
        print(f"  {lab:10s} info ratio vs L1: {r:.3f}   median |log bias|: "
              f"{s[f'{lab}_logbias'].abs().median():.3f}")
    print("\npaired within-replicate log-ratio of k-hat against L1 (0 = no change):")
    pr = []
    for lab in ("L2", "L3_greedy", "L3_nj", "L3_upgma", "L3_ens"):
        c = f"k_{lab}"
        if c not in d:
            continue
        lr = np.log(d[c].values / d["k_L1"].values)
        pr.append(dict(level=lab, median_log_ratio=float(np.median(lr)),
                       iqr=float(np.percentile(lr, 75) - np.percentile(lr, 25)),
                       frac_overestimate=float((lr > 0).mean()), n=len(lr)))
    prd = pd.DataFrame(pr); prd.to_csv(os.path.join(TAB, "e17_paired_ratios.csv"), index=False)
    print(prd.round(3).to_string(index=False))

    print("\ntopology recovery against the true sampled tree:")
    rc = [c for c in d.columns if c.startswith(("rf_recall", "rf_precision"))]
    if rc:
        print(d.groupby(["recorder", "genealogy"])[rc].median().round(3).to_string())

    print("\ntree geometry, median over replicates:")
    cols = [c for c in d.columns if c.startswith(("mrca_norm", "colless", "close_frac"))]
    print(d.groupby(["recorder", "genealogy"])[cols].median().round(3).to_string())
    print("\nrecorder state:")
    print(d.groupby("recorder")[["saturation", "missing", "n_informative_chars"]].median().round(3).to_string())


if __name__ == "__main__":
    t0 = time.time(); main(); print(f"{time.time() - t0:.0f}s")
