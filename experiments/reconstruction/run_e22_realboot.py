"""E22 (Milestone 6, Experiment 4): dataset-specific parametric bootstrap for the occupancy test.

Milestone 4 showed the chi-square reference for the shared-occupancy test is anti-conservative (14.3% type-I
error against a nominal 5%) and supplied one calibrated threshold from a 12-clone simulation. The real datasets
have different clone counts and tree shapes, so the honest calibration is dataset-specific: simulate under the
null ON THE ACTUAL TREES, refit, and read the observed statistic against that null.

Null: all clones share the pooled occupancy, each keeps its own fitted rate.
Alternative: clone-specific occupancy, clone-specific rate.

Usage: python -I run_e22_realboot.py <proc_dir> <quinn_trees> <quinn_states> <cells_tables> <out_tables> [--quick]
"""
import os, re, sys, glob
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.lineage_inference.ctmc import FastTree
from src.lineage_inference.perclone import fit_clone, lrt_shared, variance_components

QUICK = "--quick" in sys.argv
B = 40 if QUICK else 120
MIN_CELLS, MAX_CELLS, MAX_TREES = 25, 400, 12
rng = np.random.default_rng(20261205)
CPG = np.linspace(0.05, 0.95, 13)              # coarse grids for the bootstrap inner loop
CKG = np.exp(np.linspace(np.log(5e-3), np.log(20), 16))


def simulate_states(par, pi, k, leaf_mask, rng):
    """Simulate a two-state CTMC with unit edges down a fixed tree; return states at the leaves."""
    n = len(par)
    order = np.argsort([_depth(par, v) for v in range(n)])
    E = np.exp(-k)
    st = np.empty(n, int)
    for v in order:
        if par[v] < 0:
            st[v] = int(rng.random() < pi)
        else:
            p = st[par[v]]
            up = pi * (1 - E) if p == 0 else pi + (1 - pi) * E
            st[v] = int(rng.random() < up)
    out = np.full(n, -1)
    out[leaf_mask] = st[leaf_mask]
    return out


def _depth(par, v):
    c = 0
    while par[v] >= 0:
        v = par[v]; c += 1
    return c


def boot_pvalue(trees, fits, obs_stat, rng, B=B):
    """Parametric-bootstrap null for the shared-occupancy LRT, simulated on these very trees."""
    pooled = float(np.mean([f["pi"] for f in fits]))
    null = []
    for b in range(B):
        sim = []
        for t, f in zip(trees, fits):
            mask = t.leaf_state >= 0
            ls = simulate_states(t.parent, pooled, f["k"], mask, rng)
            sim.append(FastTree(t.parent, ls))
        try:
            sf = [fit_clone(t, pi_grid=CPG, k_grid=CKG, profile=False) for t in sim]
            st = lrt_shared(sim, sf, "pi", maxiter_mult=12)["stat"]
            null.append(st)
        except Exception:
            continue
    null = np.array(null)
    if len(null) < 20:
        return np.nan, np.nan, len(null)
    p = float((null >= obs_stat).mean())
    return p, float(np.percentile(null, 95)), len(null)


def simeonov_trees(proc, ctab, m):
    d = pd.read_csv(os.path.join(ctab, f"m2_real_M{m}_cells.csv.gz"))
    el = pd.read_csv(os.path.join(proc, f"M{m}/output-files/pdac_mouse{m}_full_edgelist.txt"), sep="\t")
    el["V1"] = el.V1.str.strip(); el["V2"] = el.V2.str.strip()
    state = dict(zip(d.leaf.astype(str), d.state2))
    trees, names = [], []
    for c, g in el.groupby("cloneID"):
        par = dict(zip(g.V2, g.V1))
        nodes = sorted(set(g.V1) | set(g.V2))
        obs = [n for n in nodes if n in state]
        if not (MIN_CELLS <= len(obs) <= MAX_CELLS):
            continue
        keep = set(obs)
        for n in obs:
            x = n
            while x in par:
                x = par[x]; keep.add(x)
        keep = sorted(keep); idx = {n: i for i, n in enumerate(keep)}
        parent = [idx[par[n]] if (n in par and par[n] in idx) else -1 for n in keep]
        if sum(1 for v in parent if v < 0) != 1:
            continue
        ls = [state.get(n, -1) for n in keep]
        trees.append(FastTree(parent, ls)); names.append(int(c))
    return trees, names


def quinn_trees(tdir, sdir, mouse, gsm_dir, primary_only=False):
    from experiments.quinn_analysis.run_quinn import PRIMARY, parse_newick
    f = os.path.join(sdir, f"quinn_cells_{mouse}.csv.gz")
    if not os.path.exists(f):
        return [], []
    d = pd.read_csv(f)
    d["site"] = d.sampleID.astype(str).str.replace("^30k", "", regex=True)
    if primary_only:
        d = d[d.site.isin(PRIMARY)]
        if len(d) < 200:
            return [], []
        X = np.column_stack([np.log(np.maximum(d.lib.values, 1)), d.n_genes.values, np.ones(len(d))])
        b, *_ = np.linalg.lstsq(X, d.cc_score.values, rcond=None)
        adj = d.cc_score.values - X @ b
        d = d.assign(state_use=(adj > np.median(adj)).astype(int))
    else:
        d = d.assign(state_use=d.state2_cc)
    trees, names = [], []
    for tf in sorted(glob.glob(os.path.join(tdir, gsm_dir, "trees", "*hybrid_priors*.txt"))):
        mm = re.search(r"lg(\d+)", os.path.basename(tf))
        if mm is None:
            continue
        lg = int(mm.group(1)); dl = d[d.LineageGroup == lg]
        if len(dl) < MIN_CELLS:
            continue
        s1 = dict(zip(dl.cellBC.astype(str), dl.state_use))
        s2 = dict(zip(dl.site.astype(str) + "." + dl.cellBC.astype(str), dl.state_use))
        par, nm = parse_newick(open(tf).read())
        best = None
        for sd in (s1, s2):
            ls = np.full(len(par), -1); n = 0
            for i, name in nm.items():
                if name in sd:
                    ls[i] = sd[name]; n += 1
            if best is None or n > best[1]:
                best = (ls, n)
        ls, n = best
        if MIN_CELLS <= n <= MAX_CELLS:
            trees.append(FastTree(par, ls)); names.append(lg)
    return trees, names


def run_set(label, trees, names, out):
    if len(trees) < 3:
        print(f"{label}: only {len(trees)} trees, skipped", flush=True)
        return None
    if len(trees) > MAX_TREES:
        sz = [int((t.leaf_state >= 0).sum()) for t in trees]
        keep = np.argsort(-np.array(sz))[:MAX_TREES]
        trees = [trees[i] for i in keep]; names = [names[i] for i in keep]
    fits = [fit_clone(t) for t in trees]
    vc = variance_components(fits, min_n=MIN_CELLS)
    obs = lrt_shared(trees, fits, "pi", maxiter_mult=25)
    p_boot, thr, nb = boot_pvalue(trees, fits, obs["stat"], rng)
    row = dict(dataset=label, n_clones=len(trees), lrt_stat=obs["stat"], p_chi2=obs["p_asym"],
               p_bootstrap=p_boot, boot_threshold_95=thr, n_boot=nb,
               V_pi=vc.get("V_pi", np.nan), V_pi_raw=vc.get("V_pi_raw", np.nan))
    print(f"{label}: {len(trees)} clones, LRT {obs['stat']:.1f}, chi2 p {obs['p_asym']:.2e}, "
          f"bootstrap p {p_boot:.4f} (threshold {thr:.1f}, {nb} draws), V_pi {row['V_pi']:.4f}", flush=True)
    out.append(row)
    return row


def main(proc, qtrees, qstates, ctab, out_tab):
    out = []
    for m in (1, 2):
        t, n = simeonov_trees(proc, ctab, m)
        run_set(f"Simeonov mouse{m}", t, n, out)
    GSM = {"5k": "GSM4905334_trees_5k", "30k": "GSM4905336_trees_30k",
           "10k": "GSM4905338_trees_10k", "100k": "GSM4905340_trees_100k"}
    for mouse, gd in GSM.items():
        t, n = quinn_trees(qtrees, qstates, mouse, gd, primary_only=False)
        run_set(f"Quinn {mouse} all sites", t, n, out)
        t, n = quinn_trees(qtrees, qstates, mouse, gd, primary_only=True)
        run_set(f"Quinn {mouse} primary only", t, n, out)
    d = pd.DataFrame(out)
    d.to_csv(os.path.join(out_tab, "m6_real_bootstrap_tests.csv"), index=False)
    pd.set_option("display.width", 220)
    print("\n" + d.round(4).to_string(index=False))
    sig = d[d.p_bootstrap <= 0.05]
    print(f"\neligible datasets: {len(d)}; occupancy heterogeneity significant by dataset-specific "
          f"bootstrap in {len(sig)}")
    print("done")


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    main(*a[:5])
