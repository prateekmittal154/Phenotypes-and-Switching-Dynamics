"""E19: alternative reconstructions for Simeonov from the authors' own allele-level records.

The Quinn deposit ships three published reconstructions, which E18 uses as an ensemble. Simeonov ships only one
tree, but it does ship `pdac_mouseN_final_editing_data.txt`, the per-cell per-target indel calls the tree was
built from. That lets us build independent reconstructions of the same clones and ask whether the Milestone 4
conclusions depend on the authors' particular reconstruction.

Usage: python -I run_e19_simeonov.py <proc_dir> <cells_tables_dir> <out_tables_dir>
"""
import os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.branching.recorder import reconstruct, mrca_depths, colless, MISSING
from src.lineage_inference.ctmc import VarTree, fit_var

MAX_CELLS = 220          # per clone, for tractable O(n^2)-O(n^3) reconstruction
MIN_CELLS = 30
rng = np.random.default_rng(20261204)


def allele_matrix(ed, cells):
    """cells x (intID, target) integer allele matrix; -1 for unobserved."""
    sub = ed[ed.cellID.isin(cells)]
    if len(sub) == 0:
        return None, None
    sub = sub.assign(site=sub.intID.astype(str) + "|" + sub.target.astype(str))
    codes, uniq = pd.factorize(sub.target_edit.astype(str))
    sub = sub.assign(code=codes + 1)
    piv = sub.pivot_table(index="cellID", columns="site", values="code", aggfunc="first")
    # keep sites seen in at least half the cells, then cells with enough observed sites
    piv = piv.loc[:, piv.notna().mean() >= 0.5]
    if piv.shape[1] < 3:
        return None, None
    piv = piv[piv.notna().mean(axis=1) >= 0.5]
    if len(piv) < MIN_CELLS:
        return None, None
    A = piv.fillna(MISSING).astype(int).values
    return A, list(piv.index)


def main(proc, ctab, out):
    rows, geo = [], []
    for m in (1, 2):
        ed_path = os.path.join(proc, f"M{m}/output-files/pdac_mouse{m}_final_editing_data.txt")
        cells_path = os.path.join(ctab, f"m2_real_M{m}_cells.csv.gz")
        if not (os.path.exists(ed_path) and os.path.exists(cells_path)):
            print("missing inputs for mouse", m); continue
        cells = pd.read_csv(cells_path)
        state = dict(zip(cells.leaf.astype(str), cells.state2))
        ed = pd.read_csv(ed_path, sep="\t", usecols=["cellID", "intID", "target", "target_edit", "cloneID"])
        ed["cellID"] = ed.cellID.astype(str)
        print(f"mouse {m}: editing rows {len(ed)}, clones {ed.cloneID.nunique()}", flush=True)
        for clone, g in ed.groupby("cloneID"):
            avail = [c for c in g.cellID.unique() if c in state]
            if len(avail) < MIN_CELLS:
                continue
            if len(avail) > MAX_CELLS:
                avail = list(rng.choice(avail, MAX_CELLS, replace=False))
            A, order = allele_matrix(g, set(avail))
            if A is None:
                continue
            st = np.array([state[c] for c in order])
            if st.min() == st.max():
                continue
            for method in ("greedy", "nj", "upgma"):
                par = reconstruct(A, method, rng)
                ls = np.full(len(par), -1); ls[:len(order)] = st
                t = VarTree(par, ls)
                f = fit_var(t)
                d, paths = mrca_depths(par, len(order), max_pairs=6000, rng=rng)
                h = max(len(q) for q in paths) - 1
                rows.append(dict(mouse=m, clone=int(clone), method=method, n_cells=len(order),
                                 n_sites=A.shape[1], pi=f["pi"], k=f["k"], k_at_bound=f["k_at_bound"]))
                geo.append(dict(mouse=m, clone=int(clone), method=method, n_cells=len(order), height=h,
                                mrca_norm=float(d.mean() / max(h, 1)),
                                colless=colless(par, len(order)),
                                close_frac=float((d >= 0.7 * h).mean()),
                                L_min=float(np.percentile(h - d, 10)),
                                L_max=float(np.percentile(h - d, 90))))
            print(f"  clone {clone}: {len(order)} cells, {A.shape[1]} sites", flush=True)
    r = pd.DataFrame(rows); g = pd.DataFrame(geo)
    r.to_csv(os.path.join(out, "m5_simeonov_alt_reconstruction_fits.csv"), index=False)
    g.to_csv(os.path.join(out, "m5_simeonov_alt_reconstruction_geometry.csv"), index=False)
    pd.set_option("display.width", 240)
    if len(r):
        print("\n== estimates by reconstruction method")
        print(r.groupby(["mouse", "method"]).agg(clones=("clone", "size"), pi=("pi", "median"),
                                                 k=("k", "median"),
                                                 at_bound=("k_at_bound", "mean")).round(3).to_string())
        vc = []
        for (mm, me), gg in r.groupby(["mouse", "method"]):
            if len(gg) >= 3:
                vc.append(dict(mouse=mm, method=me, n_clones=len(gg),
                               V_pi_raw=float(np.var(gg.pi, ddof=1))))
        if vc:
            print("\n== occupancy variance by reconstruction")
            print(pd.DataFrame(vc).round(4).to_string(index=False))
            pd.DataFrame(vc).to_csv(os.path.join(out, "m5_simeonov_variance_by_reconstruction.csv"),
                                    index=False)
    if len(g):
        print("\n== geometry by reconstruction")
        print(g.groupby(["mouse", "method"])[["height", "mrca_norm", "close_frac", "L_min", "L_max"]]
              .median().round(3).to_string())
    print("done")


if __name__ == "__main__":
    main(*sys.argv[1:4])
