"""E21 (Milestone 6, Experiment 3): sampling-design optimisation at fixed budget.

Five strategies, all spending the same 96 cells, across two growth models and two recorders:

  random        uniform from the final live population — what is usually done
  balanced      equal numbers per founder clone
  twotime       half at an earlier generation, half at the end
  recorder_bc   PROSPECTIVE and practical: pre-screen candidate cells by barcode similarity and preferentially
                take cells that look closely related. Uses only T2 information, so an experimenter could
                actually do it (e.g. a cheap targeted pre-screen before full transcriptome sequencing)
  oracle        picks cells with recent common ancestors using the TRUE genealogy. Unattainable; included
                only as an upper benchmark

Outcome measures: realised close-relative pair fraction, estimator bias, and profile-interval coverage for the
switching rate under reconstruction-aware inference (ensemble + ultrametric branch times, the best variant
from E20).

Usage: python -I run_e21_sampling.py [--quick]
"""
import os, sys, time
from multiprocessing import Pool
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.branching.recorder import (simulate_genealogy_tracked, evolve_states, induced_subtree,
                                    record_barcodes, reconstruct, mrca_depths, MISSING)
from src.lineage_inference.ctmc import VarTree
from experiments.reconstruction.run_e20 import ultrametric_edges, profile_k

TAB = os.path.join(ROOT, "results/tables"); os.makedirs(TAB, exist_ok=True)
QUICK = "--quick" in sys.argv
REPS = 5 if QUICK else 30
BUDGET = 96
N_GEN = 10
RECORDERS = {"gestalt_like": dict(n_sites=9, mu=0.30, n_states=12, dropout=0.01, missing=0.12),
             "cassiopeia_like": dict(n_sites=25, mu=0.12, n_states=40, dropout=0.005, missing=0.20)}
GROWTH = {"exponential": dict(p_death=0.02, bottleneck_at=None, bottleneck_size=None),
          "bottlenecked": dict(p_death=0.02, bottleneck_at=4, bottleneck_size=4)}
K_GRID = (0.3, 1.0)


def founder_of(par, v, snapshots, gen_idx=2):
    """Label each cell by its ancestor alive at an early generation, as a proxy for clone identity."""
    anc = set(snapshots[gen_idx])
    x = v
    while x >= 0 and x not in anc:
        x = par[x]
    return x


def pick(strategy, par, live, snapshots, A_all, rng):
    """Return the sampled cell indices for each strategy, all of size BUDGET."""
    live = np.asarray(live)
    if strategy == "random":
        return rng.choice(live, BUDGET, replace=False)
    if strategy == "balanced":
        fo = np.array([founder_of(par, v, snapshots) for v in live])
        groups = [live[fo == f] for f in np.unique(fo)]
        per = max(1, BUDGET // len(groups))
        out = []
        for g in groups:
            out += list(rng.choice(g, min(per, len(g)), replace=False))
        pool = np.setdiff1d(live, np.array(out))
        if len(out) < BUDGET and len(pool):
            out += list(rng.choice(pool, min(BUDGET - len(out), len(pool)), replace=False))
        return np.array(out[:BUDGET])
    if strategy == "twotime":
        early = np.asarray(snapshots[max(N_GEN - 3, 1)])
        n_e = min(BUDGET // 2, len(early))
        a = rng.choice(early, n_e, replace=False)
        b = rng.choice(live, BUDGET - n_e, replace=False)
        return np.concatenate([a, b])
    if strategy == "recorder_bc":
        # prospective: screen 4x the budget, then greedily take barcode-similar neighbours
        cand = rng.choice(live, min(4 * BUDGET, len(live)), replace=False)
        Ac = A_all[cand]
        chosen, used = [], np.zeros(len(cand), bool)
        while len(chosen) < BUDGET and not used.all():
            seeds = np.where(~used)[0]
            i = seeds[rng.integers(len(seeds))]
            m = (Ac != MISSING) & (Ac[i] != MISSING)
            sim = np.array([(Ac[j][m[j]] == Ac[i][m[j]]).mean() if m[j].sum() else 0.0
                            for j in range(len(cand))])
            sim[used] = -1
            take = np.argsort(-sim)[:min(8, BUDGET - len(chosen))]
            for t in take:
                if not used[t]:
                    used[t] = True; chosen.append(cand[t])
        return np.array(chosen[:BUDGET])
    if strategy == "oracle":
        # unattainable: pick sibling-rich sets using the true genealogy
        par_of = {v: par[v] for v in live}
        byparent = {}
        for v in live:
            byparent.setdefault(par_of[v], []).append(v)
        groups = sorted(byparent.values(), key=lambda g: -len(g))
        out = []
        for g in groups:
            out += list(g)
            if len(out) >= BUDGET:
                break
        if len(out) < BUDGET:
            pool = np.setdiff1d(live, np.array(out))
            out += list(rng.choice(pool, BUDGET - len(out), replace=False))
        return np.array(out[:BUDGET])
    raise ValueError(strategy)


def one(args):
    rec, growth, k_true, rep = args
    rng = np.random.default_rng(abs(hash((rec, growth, round(k_true, 3), rep, "e21"))) % (2 ** 32))
    rp = RECORDERS[rec]
    par, gen, live, snaps = simulate_genealogy_tracked(N_GEN, rng, **GROWTH[growth])
    if len(live) < BUDGET:
        return None
    st = evolve_states(par, rng, 0.4, k_true)
    A_all = record_barcodes(par, gen, rng, **rp)
    rows = []
    for strat in ("random", "balanced", "twotime", "recorder_bc", "oracle"):
        try:
            samp = pick(strat, par, live, snaps, A_all, rng)
        except Exception:
            continue
        samp = np.unique(samp)
        if len(samp) < 24:
            continue
        p1, l1, leafmap, _ = induced_subtree(par, samp)
        d1, paths = mrca_depths(p1, len(leafmap), max_pairs=8000, rng=rng)
        h1 = max(len(q) for q in paths) - 1
        close = float((d1 >= 0.7 * h1).mean()) if h1 > 0 else np.nan
        Al = A_all[samp]
        obs = st[samp]
        ens = []
        for m in ("greedy", "nj", "upgma"):
            p3 = reconstruct(Al, m, rng)
            ls3 = np.full(len(p3), -1); ls3[:len(samp)] = obs
            ens.append(VarTree(p3, ls3, elen=ultrametric_edges(p3, len(samp), N_GEN)))

        def ens_ll(p, k):
            v = np.array([t.loglik(p, p, k) for t in ens])
            mx = v.max()
            return float(mx + np.log(np.mean(np.exp(v - mx))))

        kh, lo, hi, _, _ = profile_k(ens_ll)
        rows.append(dict(recorder=rec, growth=growth, k_true=k_true, rep=rep, strategy=strat,
                         n_sampled=len(samp), close_frac=close, mrca_mean=float(d1.mean()),
                         k_hat=kh, lo=lo, hi=hi, cover=int(lo <= k_true <= hi),
                         log_bias=float(np.log(kh / k_true)), log_width=float(np.log(hi / lo))))
    return rows


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def main():
    jobs = [(r, g, k, rep) for r in RECORDERS for g in GROWTH for k in K_GRID for rep in range(REPS)]
    with Pool(10) as p:
        out = [r for rs in p.map(one, jobs) if rs for r in rs]
    d = pd.DataFrame(out)
    d.to_csv(os.path.join(TAB, "e21_sampling_designs.csv"), index=False)
    pd.set_option("display.width", 240)
    g = d.groupby("strategy").agg(n=("cover", "size"), close_frac=("close_frac", "median"),
                                  mrca=("mrca_mean", "median"), log_bias=("log_bias", "median"),
                                  coverage=("cover", "mean"), log_width=("log_width", "median"))
    g["cov_lo"] = [wilson(int(c * n), int(n))[0] for c, n in zip(g.coverage, g.n)]
    g["cov_hi"] = [wilson(int(c * n), int(n))[1] for c, n in zip(g.coverage, g.n)]
    g.to_csv(os.path.join(TAB, "e21_sampling_summary.csv"))
    print(g.round(3).to_string())
    print("\nby recorder and growth model, close-relative pair fraction:")
    print(d.pivot_table(index=["recorder", "growth"], columns="strategy",
                        values="close_frac").round(4).to_string())
    print("\nby recorder and growth model, interval coverage:")
    print(d.pivot_table(index=["recorder", "growth"], columns="strategy",
                        values="cover").round(3).to_string())


if __name__ == "__main__":
    t0 = time.time(); main(); print(f"{time.time() - t0:.0f}s")
