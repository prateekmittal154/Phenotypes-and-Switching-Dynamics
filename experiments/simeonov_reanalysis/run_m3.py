"""Milestone 3 re-analysis of Simeonov: per-clone occupancy and switching rate, held-out-site prediction ladder,
leave-one-clone-out, clone-1 contribution, and independent-gene validation.

Usage: python -I run_m3.py <tables_dir> <proc_dir>
"""
import json, os, sys
import numpy as np, pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.lineage_inference.perclone import fit_clone, variance_components, lrt_shared, moran_path, EXPIT, LOGIT
from src.lineage_inference.ctmc import FastTree
from src.biological_validation.extract_panel import HELDOUT
from src.model_comparison.models import build

MET = {"Met", "Liver", "Lung"}
rng = np.random.default_rng(20261130)


def load(tab, proc, m):
    d = pd.read_csv(os.path.join(tab, f"m2_real_M{m}_cells.csv.gz"))
    el = pd.read_csv(os.path.join(proc, f"M{m}/output-files/pdac_mouse{m}_full_edgelist.txt"), sep="\t")
    el["V1"] = el.V1.str.strip(); el["V2"] = el.V2.str.strip()
    return d, el


def clone_trees(d, el, sites=None, min_cells=25, drop_clone=None):
    """One FastTree per clone, restricted to the given sites."""
    sub = d if sites is None else d[d.site.isin(sites)]
    state = dict(zip(sub.leaf, sub.state2)); site = dict(zip(sub.leaf, sub.site))
    trees, names = [], []
    for c, g in el.groupby("cloneID"):
        if drop_clone is not None and c == drop_clone:
            continue
        par = dict(zip(g.V2, g.V1))
        nodes = sorted(set(g.V1) | set(g.V2))
        obs = [n for n in nodes if n in state]
        if len(obs) < min_cells:
            continue
        keep = set(obs)
        for n in obs:
            x = n
            while x in par:
                x = par[x]; keep.add(x)
        keep = sorted(keep)
        idx = {n: i for i, n in enumerate(keep)}
        parent = [idx[par[n]] if (n in par and par[n] in idx) else -1 for n in keep]
        if sum(1 for v in parent if v < 0) != 1:
            continue
        ls = [state.get(n, -1) for n in keep]
        ismet = [1 if (n in site and site[n] in MET) else 0 for n in keep]
        trees.append(build(parent, ls, ismet)); names.append(int(c))
    return trees, names


def main(tab, proc):
    M, rows_pc, rows_vc, rows_pred, rows_loo, rows_gene = {}, [], [], [], [], []
    for m in (1, 2):
        d, el = load(tab, proc, m)
        d = d[d.leaf.notna()].copy()
        trees, names = clone_trees(d, el)
        if not trees:
            continue
        fits = [fit_clone(t) for t in trees]
        for c, f in zip(names, fits):
            rows_pc.append(dict(mouse=m, clone=c, **f))
        vc = variance_components(fits)
        lp = lrt_shared(trees, fits, "pi"); lk = lrt_shared(trees, fits, "k")
        rows_vc.append(dict(mouse=m, scope="all_sites", lrt_pi_stat=lp["stat"], lrt_pi_p=lp["p_asym"],
                            lrt_k_stat=lk["stat"], lrt_k_p=lk["p_asym"], **vc))
        # primary tumour only (where founders are many; Proposition T4)
        tp, npm = clone_trees(d, el, sites=["PT"], min_cells=25)
        if len(tp) >= 3:
            fp = [fit_clone(t) for t in tp]
            vp = variance_components(fp)
            lpp = lrt_shared(tp, fp, "pi"); lkk = lrt_shared(tp, fp, "k")
            rows_vc.append(dict(mouse=m, scope="PT_only", lrt_pi_stat=lpp["stat"], lrt_pi_p=lpp["p_asym"],
                                lrt_k_stat=lkk["stat"], lrt_k_p=lkk["p_asym"], **vp))
            for c, f in zip(npm, fp):
                rows_pc.append(dict(mouse=m, clone=c, scope="PT", **f))
        # leave-one-clone-out on the variance components
        for drop in names:
            sub = [(t, f, c) for t, f, c in zip(trees, fits, names) if c != drop]
            if len(sub) >= 4:
                v = variance_components([s[1] for s in sub])
                rows_loo.append(dict(mouse=m, dropped_clone=drop, **v))

        # ---- held-out-site prediction ladder: predict per-clone high-EMT fraction in the peritoneal site
        tr, nm = clone_trees(d, el, sites=["PT", "Liver", "Lung"], min_cells=25)
        obs = d[d.site == "Met"].groupby("cloneID").state2.agg(["mean", "size"])
        obs = obs[obs["size"] >= 20]
        if len(tr) >= 3 and len(obs):
            ftr = [fit_clone(t) for t in tr]
            pi_tree = dict(zip(nm, [f["pi"] for f in ftr]))
            train = d[d.site.isin(["PT", "Liver", "Lung"])]
            grand = float(train.state2.mean())
            pt_mean = train[train.site == "PT"].groupby("cloneID").state2.mean().to_dict()
            cl_mean = train.groupby("cloneID").state2.mean().to_dict()
            freq = train.cloneID.value_counts(normalize=True).to_dict()
            for c, r in obs.iterrows():
                rows_pred.append(dict(
                    mouse=m, clone=int(c), n=int(r["size"]), observed=float(r["mean"]),
                    p_grand=grand,
                    p_clone_freq=grand,                       # clone frequency alone carries no state information
                    p_pt_mean=pt_mean.get(c, grand),
                    p_clone_mean=cl_mean.get(c, grand),
                    p_tree_pi=pi_tree.get(c, grand),
                    clone_train_n=int((train.cloneID == c).sum())))

        # ---- independent-gene validation: do held-out programmes track clone occupancy / switching rate?
        pc = pd.DataFrame([dict(clone=c, pi=f["pi"], k=f["k"], at_bound=f["k_at_bound"]) for c, f in zip(names, fits)])
        cov = d.groupby("cloneID").agg(n=("e", "size"), nfeat=("nFeature_RNA", "mean"),
                                       mt=("percent_mt", "mean"),
                                       met_frac=("site", lambda s: float(s.isin(MET).mean())))
        cc = pc.set_index("clone").join(cov)
        for prog in HELDOUT:
            pv = d.groupby("cloneID")[prog].mean()
            g = cc.join(pv.rename("prog")).dropna(subset=["prog"])
            if len(g) >= 5:
                r_pi = stats.spearmanr(g.pi, g.prog)
                # partial correlation controlling clone size, depth, mito and compartment mix
                X = np.column_stack([np.log(g.n), g.nfeat, g.mt, g.met_frac, np.ones(len(g))])
                def resid(y):
                    b, *_ = np.linalg.lstsq(X, np.asarray(y, float), rcond=None)
                    return np.asarray(y, float) - X @ b
                rp = stats.spearmanr(resid(g.pi), resid(g.prog))
                gk = g[~g.at_bound]
                r_k = stats.spearmanr(np.log(gk.k), gk.prog) if len(gk) >= 5 else None
                rows_gene.append(dict(mouse=m, prog=prog, n_clones=len(g),
                                      rho_pi=float(r_pi.statistic), p_pi=float(r_pi.pvalue),
                                      rho_pi_adj=float(rp.statistic), p_pi_adj=float(rp.pvalue),
                                      n_k=len(gk),
                                      rho_logk=float(r_k.statistic) if r_k else np.nan,
                                      p_logk=float(r_k.pvalue) if r_k else np.nan))
        M[f"M{m}_n_clone_trees"] = len(trees)
        print("mouse", m, "done", flush=True)

    pd.DataFrame(rows_pc).to_csv(os.path.join(tab, "m3_real_perclone.csv"), index=False)
    pd.DataFrame(rows_vc).to_csv(os.path.join(tab, "m3_real_variance_components.csv"), index=False)
    pd.DataFrame(rows_loo).to_csv(os.path.join(tab, "m3_real_leave_one_clone_out.csv"), index=False)
    pd.DataFrame(rows_gene).to_csv(os.path.join(tab, "m3_real_gene_validation.csv"), index=False)
    p = pd.DataFrame(rows_pred)
    if len(p):
        for col in ("p_grand", "p_pt_mean", "p_clone_mean", "p_tree_pi"):
            p["err_" + col] = (p[col] - p.observed).abs()
        p.to_csv(os.path.join(tab, "m3_real_heldout_ladder.csv"), index=False)
        print(p.groupby("mouse")[[c for c in p.columns if c.startswith("err_")]].mean().round(4).to_string())
        # excluding the dominant clone 1 / clone 2
        for m_ in p.mouse.unique():
            q = p[(p.mouse == m_) & (p.n < p[p.mouse == m_].n.max())]
            if len(q):
                print(f"  mouse {m_} excluding its largest clone (n={len(q)} clones):",
                      q[[c for c in p.columns if c.startswith("err_")]].mean().round(4).to_dict())
    with open(os.path.join(tab, "..", "metrics", "m3_real_metrics.json"), "w") as fh:
        json.dump(M, fh, indent=1, default=float)
    print("done")


if __name__ == "__main__":
    main(*sys.argv[1:3])
