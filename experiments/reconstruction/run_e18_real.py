"""E18: reconstruction uncertainty in the real data.

The Quinn deposit contains THREE published reconstructions of every lineage group — `hybrid_priors.alleleThresh`
(the authors' primary), `nj` and `nj_priors`. These are defensible alternative trees produced by the original
authors from the same alleles, so they give a genuine reconstruction ensemble at no modelling cost.

For every clone we refit (pi_c, k_c) on each reconstruction, compare geometry and estimates across them, and
form an ensemble estimate by averaging the likelihood over trees. The question is whether the Milestone 4
conclusions - clone occupancy heterogeneous, switching rates unidentifiable - survive reconstruction
uncertainty.

Usage: python -I run_e18_real.py <quinn_dir> <states_dir> <tables_dir>
"""
import os, re, sys, glob
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.lineage_inference.ctmc import VarTree, fit_var
from src.lineage_inference.perclone import variance_components
from src.branching.recorder import mrca_depths, colless
from experiments.quinn_analysis.run_quinn import TREES, PRIMARY, parse_newick

MAX_LEAVES = 1200
MAX_TREES = 25
VARIANTS = {"hybrid": "*hybrid_priors*", "nj": "*_tree_nj.processed*", "nj_priors": "*_tree_nj_priors*"}
rng = np.random.default_rng(20261203)


def build_tree(nwk, state_of):
    par, names = parse_newick(nwk)
    ls = np.full(len(par), -1)
    n = 0
    leaf_idx = []
    for idx, nm in names.items():
        if nm in state_of:
            ls[idx] = state_of[nm]; n += 1; leaf_idx.append(idx)
    if n < 25 or n > MAX_LEAVES:
        return None, n
    return VarTree(par, ls), n


def main(qdir, sdir, tab):
    rows, geo = [], []
    for mouse, tdir in TREES.items():
        f = os.path.join(sdir, f"quinn_cells_{mouse}.csv.gz")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        d["site"] = d.sampleID.astype(str).str.replace("^30k", "", regex=True)
        per_variant = {}
        for vname, pat in VARIANTS.items():
            files = sorted(glob.glob(os.path.join(qdir, tdir, "trees", pat + ".txt")))
            got = {}
            for tf in files:
                m = re.search(r"lg(\d+)", os.path.basename(tf))
                if m is None:
                    continue
                lg = int(m.group(1)); dl = d[d.LineageGroup == lg]
                if len(dl) < 25:
                    continue
                s1 = dict(zip(dl.cellBC.astype(str), dl.state2_cc))
                s2 = dict(zip(dl.site.astype(str) + "." + dl.cellBC.astype(str), dl.state2_cc))
                nwk = open(tf).read()
                t1, n1 = build_tree(nwk, s1); t2, n2 = build_tree(nwk, s2)
                t, n = (t1, n1) if n1 >= n2 else (t2, n2)
                if t is not None:
                    got[lg] = (t, n)
            per_variant[vname] = got
            print(f"{mouse} {vname}: {len(got)} usable trees", flush=True)
        common = set.intersection(*[set(v) for v in per_variant.values()]) if per_variant else set()
        common = sorted(common, key=lambda c: -per_variant["hybrid"][c][1])[:MAX_TREES]
        print(f"{mouse}: {len(common)} clones present in all three reconstructions", flush=True)
        fits = {v: [] for v in VARIANTS}
        for lg in common:
            trees = {v: per_variant[v][lg][0] for v in VARIANTS}
            n = per_variant["hybrid"][lg][1]
            for v, t in trees.items():
                fr = fit_var(t)
                fits[v].append(dict(pi=fr["pi"], k=fr["k"], n=n, k_at_bound=fr["k_at_bound"],
                                    pi_se=np.nan, logk_se=np.nan))
                nl = int((t.leaf_state >= 0).sum())
                dd, paths = mrca_depths(t.parent, len(t.parent), max_pairs=4000, rng=rng)
                rows.append(dict(mouse=mouse, clone=lg, variant=v, n_leaves=nl,
                                 pi=fr["pi"], k=fr["k"], k_at_bound=fr["k_at_bound"]))
            # ensemble: average the likelihood across the three published reconstructions
            best = (np.inf, np.nan, np.nan)
            for p in np.linspace(0.08, 0.92, 15):
                for k in np.exp(np.linspace(np.log(5e-3), np.log(15), 20)):
                    v_ = [t.loglik(p, p, k) for t in trees.values()]
                    mx = max(v_)
                    nll = -(mx + np.log(np.mean(np.exp(np.array(v_) - mx))))
                    if nll < best[0]:
                        best = (nll, p, k)
            rows.append(dict(mouse=mouse, clone=lg, variant="ensemble", n_leaves=n,
                             pi=best[1], k=best[2], k_at_bound=bool(best[2] <= 6e-3 or best[2] >= 14.9)))
        # geometry per variant
        for v in VARIANTS:
            for lg in common:
                t, n = per_variant[v][lg]
                nl = int((t.leaf_state >= 0).sum())
                leaves = np.where(t.leaf_state >= 0)[0]
                par = t.parent
                paths = []
                for x in leaves:
                    p, y = [], x
                    while y >= 0:
                        p.append(y); y = par[y]
                    paths.append(tuple(reversed(p)))
                h = max(len(q) for q in paths) - 1
                pick = rng.choice(len(paths), min(len(paths), 90), replace=False)
                dd = []
                for a in range(len(pick)):
                    for b in range(a + 1, len(pick)):
                        x, y = paths[pick[a]], paths[pick[b]]
                        s = 0
                        for u, w in zip(x, y):
                            if u != w:
                                break
                            s += 1
                        dd.append(s)
                dd = np.array(dd)
                geo.append(dict(mouse=mouse, clone=lg, variant=v, n_leaves=nl, height=h,
                                mrca_mean=float(dd.mean()), mrca_norm=float(dd.mean() / max(h, 1)),
                                close_frac=float((dd >= 0.7 * h).mean()),
                                L_min=float(np.percentile(h - dd, 10)),
                                L_max=float(np.percentile(h - dd, 90))))
    r = pd.DataFrame(rows); g = pd.DataFrame(geo)
    r.to_csv(os.path.join(tab, "m5_quinn_reconstruction_fits.csv"), index=False)
    g.to_csv(os.path.join(tab, "m5_quinn_reconstruction_geometry.csv"), index=False)
    if len(r):
        pd.set_option("display.width", 240)
        print("\n== per-clone estimates by reconstruction (median over clones)")
        print(r.groupby(["mouse", "variant"]).agg(clones=("clone", "size"), pi=("pi", "median"),
                                                  k=("k", "median"),
                                                  at_bound=("k_at_bound", "mean")).round(3).to_string())
        print("\n== occupancy variance by reconstruction")
        vc = []
        for (m, v), gg in r.groupby(["mouse", "variant"]):
            if len(gg) >= 3:
                vc.append(dict(mouse=m, variant=v, n_clones=len(gg),
                               V_pi_raw=float(np.var(gg.pi, ddof=1)),
                               V_logk_raw=float(np.var(np.log(gg[~gg.k_at_bound].k), ddof=1))
                               if (~gg.k_at_bound).sum() >= 3 else np.nan))
        v = pd.DataFrame(vc); v.to_csv(os.path.join(tab, "m5_quinn_variance_by_reconstruction.csv"), index=False)
        print(v.round(4).to_string(index=False))
    if len(g):
        print("\n== geometry by reconstruction (median)")
        print(g.groupby(["mouse", "variant"])[["height", "mrca_norm", "close_frac", "L_min", "L_max"]]
              .median().round(3).to_string())
    print("done")


if __name__ == "__main__":
    main(*sys.argv[1:4])
