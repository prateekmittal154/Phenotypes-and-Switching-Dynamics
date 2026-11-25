"""E20 (Milestone 6, Experiments 1 and 2): reconstruction-aware interval calibration and branch-time models.

Six inference variants on identical observed leaf states, so each comparison isolates one error source:

  V1  true sampled topology, TRUE branch times                      reference
  V2  true topology, unit edges                                     branch-time miscalibration only
  V3  reconstructed topology, unit edges                            + topology error (the usual practice)
  V4  reconstructed topology, SCALED edges                          + a global time calibration
  V5  ensemble over three reconstructions, scaled edges             + reconstruction uncertainty
  V6  as V5 but integrating over an uncertain time scale            + branch-time uncertainty

The scale correction in V4-V6 rests on an identity worth stating: multiplying every branch length by c and
dividing k by c leaves the likelihood unchanged, so with unit edges the rate is identified only up to the
unknown number of divisions per resolved branch. If the experiment knows roughly how many generations D
elapsed, and the reconstructed tree has height h, then c = D / h is an estimate of that factor and is
information the experiment already has. V6 treats c as uncertain (log-normal, sd 0.3) and integrates.

Profile-likelihood intervals are computed for every variant and their coverage measured, at rates below,
within and above the empirically identifiable window.

Usage: python -I run_e20.py [--quick]
"""
import os, sys, time
from multiprocessing import Pool
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.branching.recorder import (simulate_genealogy, evolve_states, induced_subtree, record_barcodes,
                                    reconstruct, MISSING)
from src.lineage_inference.ctmc import VarTree

TAB = os.path.join(ROOT, "results/tables"); os.makedirs(TAB, exist_ok=True)
QUICK = "--quick" in sys.argv
REPS = 8 if QUICK else 60
N_SAMPLE = 96
N_GEN = 10
# below / within / above the empirically identifiable window for these trees
K_GRID = (0.05, 0.3, 1.0, 3.0)
RECORDERS = {"gestalt_like": dict(n_sites=9, mu=0.30, n_states=12, dropout=0.01, missing=0.12),
             "cassiopeia_like": dict(n_sites=25, mu=0.12, n_states=40, dropout=0.005, missing=0.20)}
KG = np.exp(np.linspace(np.log(3e-3), np.log(25), 36))
PG = np.linspace(0.04, 0.96, 20)
CHI1 = 3.841458820694124      # chi2(1) at 0.95, for profile intervals


def ultrametric_edges(par, n_leaves, total_time):
    """Branch times making every root-to-leaf path equal `total_time`.

    A single global scale c = D / h is NOT a valid calibration: it only makes the DEEPEST path equal D, so the
    many shallower leaves get paths that are too short and the rate is inflated. (E20 measures this failure.)
    Setting node times by t(v) = D * depth(v) / (depth(v) + max_leaf_depth_below(v)) makes the tree ultrametric,
    so every leaf sits at time D, which is what a real experiment knows: all sampled cells are the same age.
    """
    par = np.asarray(par); n = len(par)
    ch = {}
    for v, pa in enumerate(par):
        if pa >= 0:
            ch.setdefault(pa, []).append(v)
    depth = np.zeros(n, int)
    order = np.argsort([_d(par, v) for v in range(n)])
    for v in order:
        depth[v] = 0 if par[v] < 0 else depth[par[v]] + 1
    below = np.zeros(n, int)
    for v in order[::-1]:
        below[v] = 0 if v not in ch else 1 + max(below[u] for u in ch[v])
    t = np.where(depth + below > 0, total_time * depth / np.maximum(depth + below, 1), 0.0)
    el = np.zeros(n)
    for v in range(n):
        el[v] = 0.0 if par[v] < 0 else max(t[v] - t[par[v]], 1e-6)
    return el


def _d(par, v):
    c = 0
    while par[v] >= 0:
        v = par[v]; c += 1
    return c


def profile_k(loglik_fn, kgrid=KG, pgrid=PG):
    """Profile likelihood over pi for each k. Returns k-hat, the 95% profile interval and the bound flags."""
    prof = np.array([max(loglik_fn(p, k) for p in pgrid) for k in kgrid])
    i = int(np.argmax(prof))
    cut = prof[i] - CHI1 / 2
    ok = np.where(prof >= cut)[0]
    lo, hi = kgrid[ok.min()], kgrid[ok.max()]
    return float(kgrid[i]), float(lo), float(hi), bool(ok.min() == 0), bool(ok.max() == len(kgrid) - 1)


def one(args):
    rec_name, k_true, rep = args
    rng = np.random.default_rng(abs(hash((rec_name, round(k_true, 4), rep, "m6"))) % (2 ** 32))
    rp = RECORDERS[rec_name]
    pi_true = 0.4
    par, gen, live = simulate_genealogy(N_GEN, rng, p_death=0.02)
    if len(live) < N_SAMPLE:
        return None
    st = evolve_states(par, rng, pi_true, k_true)
    samp = rng.choice(live, N_SAMPLE, replace=False)
    p1, l1, leafmap, _ = induced_subtree(par, samp)
    ls1 = np.full(len(p1), -1)
    for ni, orig in leafmap.items():
        ls1[ni] = st[orig]
    obs = np.array([st[s] for s in samp])
    A = record_barcodes(par, gen, rng, **rp)[samp]

    out = dict(recorder=rec_name, k_true=k_true, rep=rep,
               missing=float((A == MISSING).mean()))
    # ---- V1 and V2 on the true topology
    for lab, el in (("V1", l1), ("V2", np.ones(len(p1)))):
        t = VarTree(p1, ls1, elen=el)
        kh, lo, hi, bl, bh = profile_k(lambda p, k: t.loglik(p, p, k))
        out.update({f"k_{lab}": kh, f"lo_{lab}": lo, f"hi_{lab}": hi,
                    f"cov_{lab}": int(lo <= k_true <= hi), f"w_{lab}": np.log(hi / lo)})
    # ---- reconstructions
    trees, heights = {}, {}
    for m in ("greedy", "nj", "upgma"):
        p3 = reconstruct(A, m, rng)
        ls3 = np.full(len(p3), -1); ls3[:N_SAMPLE] = obs
        # tree height in resolved branches
        h = 0
        for v in range(N_SAMPLE):
            d, x = 0, v
            while p3[x] >= 0:
                x = p3[x]; d += 1
            h = max(h, d)
        trees[m] = (p3, ls3); heights[m] = max(h, 1)
    # V3: reconstructed, unit edges (the usual practice), using NJ as the single tree
    p3, ls3 = trees["nj"]
    t3 = VarTree(p3, ls3)
    kh, lo, hi, bl, bh = profile_k(lambda p, k: t3.loglik(p, p, k))
    out.update(k_V3=kh, lo_V3=lo, hi_V3=hi, cov_V3=int(lo <= k_true <= hi), w_V3=np.log(hi / lo))
    # V4a: naive global scale c = D / h (kept because it FAILS, and that is informative)
    c_nj = N_GEN / heights["nj"]
    t4a = VarTree(p3, ls3, elen=np.full(len(p3), c_nj))
    kh, lo, hi, bl, bh = profile_k(lambda p, k: t4a.loglik(p, p, k))
    out.update(k_V4a=kh, lo_V4a=lo, hi_V4a=hi, cov_V4a=int(lo <= k_true <= hi), w_V4a=np.log(hi / lo))
    # V4: reconstructed, ULTRAMETRIC branch times (every leaf at age D)
    t4 = VarTree(p3, ls3, elen=ultrametric_edges(p3, N_SAMPLE, N_GEN))
    kh, lo, hi, bl, bh = profile_k(lambda p, k: t4.loglik(p, p, k))
    out.update(k_V4=kh, lo_V4=lo, hi_V4=hi, cov_V4=int(lo <= k_true <= hi), w_V4=np.log(hi / lo),
               scale_nj=c_nj)
    # V5: ensemble over the three reconstructions, each scaled by its own c
    ens = [VarTree(trees[m][0], trees[m][1], elen=ultrametric_edges(trees[m][0], N_SAMPLE, N_GEN))
           for m in trees]

    def ens_ll(p, k):
        v = np.array([t.loglik(p, p, k) for t in ens])
        mx = v.max()
        return float(mx + np.log(np.mean(np.exp(v - mx))))

    kh, lo, hi, bl, bh = profile_k(ens_ll)
    out.update(k_V5=kh, lo_V5=lo, hi_V5=hi, cov_V5=int(lo <= k_true <= hi), w_V5=np.log(hi / lo))
    # V6: ensemble, integrating over an uncertain time scale (log-normal, sd 0.3)
    cs = np.exp(np.array([-0.45, 0.0, 0.45]))
    wts = np.array([0.25, 0.5, 0.25])
    ens6 = {m: [VarTree(trees[m][0], trees[m][1],
                        elen=ultrametric_edges(trees[m][0], N_SAMPLE, N_GEN * c)) for c in cs]
            for m in trees}

    def ens6_ll(p, k):
        v = []
        for m in ens6:
            for w, t in zip(wts, ens6[m]):
                v.append(np.log(w) + t.loglik(p, p, k))
        v = np.array(v); mx = v.max()
        return float(mx + np.log(np.sum(np.exp(v - mx)) / len(trees)))

    kh, lo, hi, bl, bh = profile_k(ens6_ll)
    out.update(k_V6=kh, lo_V6=lo, hi_V6=hi, cov_V6=int(lo <= k_true <= hi), w_V6=np.log(hi / lo))
    return out


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def main():
    jobs = [(r, k, rep) for r in RECORDERS for k in K_GRID for rep in range(REPS)]
    with Pool(12) as p:
        rows = [x for x in p.map(one, jobs) if x]
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(TAB, "e20_interval_calibration.csv"), index=False)
    recs = []
    for (r, kt), g in d.groupby(["recorder", "k_true"]):
        for v in ("V1", "V2", "V3", "V4a", "V4", "V5", "V6"):
            n = len(g); c = int(g[f"cov_{v}"].sum())
            lo, hi = wilson(c, n)
            recs.append(dict(recorder=r, k_true=kt, variant=v, n=n,
                             coverage=c / n, cov_lo=lo, cov_hi=hi,
                             log_bias=float(np.median(np.log(g[f"k_{v}"] / kt))),
                             median_log_width=float(np.median(g[f"w_{v}"]))))
    s = pd.DataFrame(recs)
    s.to_csv(os.path.join(TAB, "e20_calibration_summary.csv"), index=False)
    pd.set_option("display.width", 240)
    print(s.pivot_table(index=["recorder", "k_true"], columns="variant",
                        values="coverage").round(3).to_string())
    print("\nmedian log bias:")
    print(s.pivot_table(index=["recorder", "k_true"], columns="variant",
                        values="log_bias").round(3).to_string())
    print("\nmedian log interval width:")
    print(s.pivot_table(index=["recorder", "k_true"], columns="variant",
                        values="median_log_width").round(3).to_string())
    print("\npooled over settings (Monte Carlo uncertainty from Wilson intervals):")
    for v in ("V1", "V2", "V3", "V4a", "V4", "V5", "V6"):
        c = int(d[f"cov_{v}"].sum()); n = len(d); lo, hi = wilson(c, n)
        print(f"  {v}: coverage {c / n:.3f} [{lo:.3f}, {hi:.3f}]  "
              f"median log bias {np.median(np.log(d[f'k_{v}'] / d.k_true)):+.3f}  "
              f"median log width {np.median(d[f'w_{v}']):.3f}")


if __name__ == "__main__":
    t0 = time.time(); main(); print(f"{time.time() - t0:.0f}s")
