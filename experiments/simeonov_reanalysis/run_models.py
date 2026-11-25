"""Model comparison M0-M4 on the real Simeonov within-clone trees, plus a held-out-site predictive check.

Usage: python -I run_models.py <tables_dir> <proc_dir>
"""
import json, os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.model_comparison.models import build, fit, EXPIT, LOGIT

MET = {"Met", "Liver", "Lung"}


def trees_for(m, tab, proc, min_cells=40, drop_sites=()):
    d = pd.read_csv(os.path.join(tab, f"m2_real_M{m}_cells.csv.gz"))
    d = d[~d.site.isin(drop_sites)]
    el = pd.read_csv(os.path.join(proc, f"M{m}/output-files/pdac_mouse{m}_full_edgelist.txt"), sep="\t")
    el["V1"] = el.V1.str.strip(); el["V2"] = el.V2.str.strip()
    state = dict(zip(d.leaf, d.state2)); site = dict(zip(d.leaf, d.site))
    out, names = [], []
    for c, sub in el.groupby("cloneID"):
        nodes = sorted(set(sub.V1) | set(sub.V2))
        obs = [n for n in nodes if n in state]
        if len(obs) < min_cells:
            continue
        keep = set(obs)
        par = dict(zip(sub.V2, sub.V1))
        for n in obs:                                   # keep ancestors of observed leaves
            x = n
            while x in par:
                x = par[x]; keep.add(x)
        keep = [n for n in nodes if n in keep] + [r for r in set(sub.V1) - set(sub.V2) if r in keep]
        keep = sorted(set(keep))
        idx = {n: i for i, n in enumerate(keep)}
        parent = [idx.get(par.get(n, None), -1) if par.get(n, None) in idx else -1 for n in keep]
        ls = [state.get(n, -1) for n in keep]
        isx = [1 if (n in site and site[n] in MET) else 0 for n in keep]
        if sum(1 for v in parent if v < 0) != 1:        # keep single-rooted trees only
            continue
        out.append(build(parent, ls, isx)); names.append(int(c))
    return out, names, d


def main(tab, proc):
    rows, pred = [], []
    for m in (1, 2):
        trees, names, d = trees_for(m, tab, proc, drop_sites=("PTab", "Blood"))
        if not trees:
            continue
        fits = {mod: fit(trees, mod) for mod in ("M0", "M1", "M2", "M4")}
        best = min(fits.values(), key=lambda f: f["aic"])
        for mod, f in fits.items():
            rows.append(dict(mouse=m, n_clone_trees=len(trees),
                             n_leaves=int(sum((t.leaf_state >= 0).sum() for t in trees)),
                             delta_aic=f["aic"] - best["aic"],
                             **{k: v for k, v in f.items() if k != "x"}))
        if "M4" in fits:
            x = np.array(fits["M4"]["x"])
            print(f"M{m} M4: k={np.exp(x[0]):.3f} met_logit_offset={x[1]:+.3f} "
                  f"clone pi range {EXPIT(x[2:]).min():.2f}-{EXPIT(x[2:]).max():.2f}", flush=True)
        # held-out site: fit without peritoneal leaves, predict each clone's peritoneal high-state fraction
        tr2, nm2, _ = trees_for(m, tab, proc, drop_sites=("PTab", "Blood", "Met"))
        if tr2:
            f0, f1 = fit(tr2, "M0"), fit(tr2, "M1")
            pi0 = EXPIT(np.array(f0["x"])[0])
            pi1 = dict(zip(nm2, EXPIT(np.array(f1["x"])[1:])))
            obs = d[d.site == "Met"].groupby("cloneID").state2.agg(["mean", "size"])
            obs = obs[obs["size"] >= 20]
            for c, r in obs.iterrows():
                if c in pi1:
                    pred.append(dict(mouse=m, clone=int(c), n=int(r["size"]), observed=float(r["mean"]),
                                     pred_M0=float(pi0), pred_M1=float(pi1[c])))
        print("mouse", m, "done", flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(tab, "m2_real_model_comparison.csv"), index=False)
    p = pd.DataFrame(pred)
    if len(p):
        for mod in ("M0", "M1"):
            p[f"err_{mod}"] = (p[f"pred_{mod}"] - p.observed).abs()
        p.to_csv(os.path.join(tab, "m2_real_heldout_site.csv"), index=False)
        print(p.groupby("mouse")[["err_M0", "err_M1"]].mean().round(4).to_string())
    print("done")


if __name__ == "__main__":
    main(*sys.argv[1:3])
