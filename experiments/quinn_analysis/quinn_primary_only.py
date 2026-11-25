"""Quinn occupancy/rate estimation restricted to PRIMARY-SITE cells only.

The technical control shows that between-site comparisons in this dataset are confounded by library complexity
(genes detected shifts up to 1.4 SD within clone, more than the cell-cycle signal itself). Restricting to the
primary site removes the between-tissue technical contrast, so if clone-associated occupancy survives here it
is not a between-site artefact. Depth is additionally regressed out of the state before binarising.

Usage: python -I quinn_primary_only.py <quinn_dir> <states_dir> <tables_dir>
"""
import os, re, sys, glob
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.lineage_inference.perclone import fit_clone, variance_components, lrt_shared
from experiments.quinn_analysis.run_quinn import TREES, PRIMARY, build, MAX_LEAVES, MAX_TREES


def main(qdir, sdir, tab):
    rows_vc, rows_pc = [], []
    for mouse, tdir in TREES.items():
        f = os.path.join(sdir, f"quinn_cells_{mouse}.csv.gz")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        d["site"] = d.sampleID.astype(str).str.replace("^30k", "", regex=True)
        d = d[d.site.isin(PRIMARY)].copy()
        if len(d) < 200:
            print("too few primary cells for", mouse); continue
        # residualise the cell-cycle score on log library size and genes detected, then binarise at the median
        X = np.column_stack([np.log(np.maximum(d.lib.values, 1)), d.n_genes.values, np.ones(len(d))])
        y = d.cc_score.values
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        d["cc_adj"] = y - X @ b
        d["state_adj"] = (d.cc_adj > np.median(d.cc_adj)).astype(int)
        trees, names, nmatch = [], [], []
        for tf in sorted(glob.glob(os.path.join(qdir, tdir, "trees", "*hybrid_priors*.txt"))):
            lg = re.search(r"lg(\d+)", os.path.basename(tf))
            if lg is None:
                continue
            lgi = int(lg.group(1)); dl = d[d.LineageGroup == lgi]
            if len(dl) < 25:
                continue
            nwk = open(tf).read()
            s1 = dict(zip(dl.cellBC.astype(str), dl.state_adj))
            s2 = dict(zip(dl.site.astype(str) + "." + dl.cellBC.astype(str), dl.state_adj))
            t1, n1 = build(nwk, s1); t2, n2 = build(nwk, s2)
            t, n = (t1, n1) if n1 >= n2 else (t2, n2)
            if t is not None:
                trees.append(t); names.append(lgi); nmatch.append(n)
        if len(trees) > MAX_TREES:
            keep = np.argsort(-np.array(nmatch))[:MAX_TREES]
            trees = [trees[i] for i in keep]; names = [names[i] for i in keep]; nmatch = [nmatch[i] for i in keep]
        print(mouse, "primary-only trees", len(trees), "cells", sum(nmatch), flush=True)
        if len(trees) < 3:
            continue
        fits = [fit_clone(t) for t in trees]
        for c, fi, nm in zip(names, fits, nmatch):
            rows_pc.append(dict(mouse=mouse, clone=c, n_matched=nm, **fi))
        vc = variance_components(fits, min_n=25)
        lp = lrt_shared(trees, fits, "pi", maxiter_mult=25)
        lk = lrt_shared(trees, fits, "k", maxiter_mult=25)
        rows_vc.append(dict(mouse=mouse, scope="primary_only_depth_adjusted",
                            lrt_pi_stat=lp["stat"], lrt_pi_p=lp["p_asym"],
                            lrt_k_stat=lk["stat"], lrt_k_p=lk["p_asym"], **vc))
    pd.DataFrame(rows_vc).to_csv(os.path.join(tab, "m4_quinn_primary_only.csv"), index=False)
    pd.DataFrame(rows_pc).to_csv(os.path.join(tab, "m4_quinn_primary_only_perclone.csv"), index=False)
    if rows_vc:
        print(pd.DataFrame(rows_vc).round(4).to_string())
    print("done")


if __name__ == "__main__":
    main(*sys.argv[1:4])
