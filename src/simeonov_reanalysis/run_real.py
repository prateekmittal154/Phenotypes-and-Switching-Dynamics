"""Simeonov (GSE173958) re-analysis, following protocols/m2_simeonov_prereg.md.

Usage: python -I run_real.py <panel_dir> <proc_dir> <out_tables> <out_metrics>
Confirmatory analyses H1-H5 plus the pre-specified exploratory items. Mice are analysed separately; clones are the
bootstrap unit.
"""
import json, os, sys
import numpy as np, pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.biological_validation.extract_panel import STATE_E, STATE_M, HELDOUT
from src.lineage_inference.propensity import decompose, curve_frame, shared_depth

SITES = {"PT": "Primary tumour", "Met": "Peritoneal met", "Liver": "Liver met", "Lung": "Lung met",
         "PTab": "Surgical site", "Blood": "CTC"}
MET_SITES = ["Met", "Liver", "Lung"]
rng_global = np.random.default_rng(20261126)


def emt_score(d):
    """Per-gene z-scores computed on the primary-tumour cells of this mouse, then mean(M) - mean(E)."""
    pt = d.site == "PT"
    z = {}
    for g in STATE_E + STATE_M:
        mu, sd = d.loc[pt, g].mean(), d.loc[pt, g].std()
        z[g] = (d[g] - mu) / sd if sd > 1e-9 else 0.0
    Z = pd.DataFrame(z, index=d.index)
    return Z[STATE_M].mean(axis=1) - Z[STATE_E].mean(axis=1)


def prog_scores(d):
    out = {}
    pt = d.site == "PT"
    for name, genes in HELDOUT.items():
        gs = [g for g in genes if g in d.columns]
        Z = pd.DataFrame({g: (d[g] - d.loc[pt, g].mean()) / (d.loc[pt, g].std() or 1) for g in gs}, index=d.index)
        out[name] = Z.mean(axis=1)
    return pd.DataFrame(out, index=d.index)


def load(panel_dir, proc_dir, m):
    d = pd.read_csv(os.path.join(panel_dir, f"panel_M{m}.csv.gz"))
    d["site"] = d.rowname.str.split("_").str[1]
    d["leaf"] = d.rowname.str.split("_", n=1).str[1]
    d["e"] = emt_score(d)
    d = pd.concat([d, prog_scores(d)], axis=1)
    e = pd.read_csv(os.path.join(proc_dir, f"M{m}/output-files/pdac_mouse{m}_full_edgelist.txt"), sep="\t")
    e["V1"] = e.V1.str.strip(); e["V2"] = e.V2.str.strip()
    par = dict(zip(e.V2, e.V1))
    paths = {}
    for leaf in set(e.V2) - set(e.V1):
        p, n = [], leaf
        while n in par:
            n = par[n]
            if n != "Cancer_root":
                p.append(n)
        paths[leaf] = tuple(reversed(p))
    d["path"] = d.leaf.map(paths)
    return d, e


def clone_boot(vals, clones, stat, n=2000, rng=None):
    """Clone-level bootstrap CI for a statistic computed on a (values, clones) pair."""
    rng = rng or rng_global
    uc = np.unique(clones); out = []
    for _ in range(n):
        pick = rng.choice(uc, len(uc), replace=True)
        idx = np.concatenate([np.where(clones == c)[0] for c in pick])
        v = stat(vals[idx], clones[idx])
        if v is not None and np.isfinite(v):
            out.append(v)
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))) if len(out) > 20 else (np.nan, np.nan)


def main(panel_dir, proc_dir, tab, met):
    os.makedirs(tab, exist_ok=True); os.makedirs(met, exist_ok=True)
    M = {}
    rows_h1, rows_h3, rows_dec, rows_curve, rows_prog, rows_h2 = [], [], [], [], [], []
    for m in (1, 2):
        d, el = load(panel_dir, proc_dir, m)
        d = d[d.path.notna()].copy()
        qs = d.loc[d.site == "PT", "e"].quantile([.25, .5, .75]).values
        d["state4"] = np.digitize(d.e, qs)                      # 0 E, 1 EH, 2 LH, 3 M
        d["state2"] = (d.e > qs[1]).astype(int)
        M[f"M{m}_cells_with_tree"] = int(len(d))
        M[f"M{m}_clones"] = int(d.cloneID.nunique())
        M[f"M{m}_pt_quartile_cuts"] = [float(x) for x in qs]
        if "pseudotime" in d.columns and d.pseudotime.notna().any():
            r = stats.spearmanr(d.e, d.pseudotime, nan_policy="omit")
            M[f"M{m}_emt_vs_authors_pseudotime_spearman"] = [float(r.statistic), float(r.pvalue)]
        if "adj_aggro" in d.columns:
            cl = d.groupby("cloneID").agg(aggro=("adj_aggro", "first"), n=("e", "size"),
                                          pt_n=("site", lambda s: int((s == "PT").sum())))
            ptlh = d[d.site == "PT"].groupby("cloneID").apply(
                lambda g: float((g.state4 == 2).mean()), include_groups=False).rename("pt_frac_LH")
            cl = cl.join(ptlh)
            sel = cl[(cl.pt_n >= 10) & cl.aggro.notna() & cl.pt_frac_LH.notna()]
            if len(sel) >= 5:
                r = stats.spearmanr(sel.pt_frac_LH, sel.aggro)
                M[f"M{m}_H1_clone_PTfracLH_vs_aggression_spearman"] = [float(r.statistic), float(r.pvalue), len(sel)]
            cl.to_csv(os.path.join(tab, f"m2_real_M{m}_clones.csv"))

        # ---- H1: site enrichment, overall and within clone
        base = d[d.site == "PT"]
        for s in MET_SITES:
            sub = d[d.site == s]
            if len(sub) < 50:
                continue
            raw = float(sub.e.mean() - base.e.mean())
            shared = sorted(set(sub.cloneID) & set(base.cloneID))
            per = [(c, float(sub[sub.cloneID == c].e.mean() - base[base.cloneID == c].e.mean()),
                    int((sub.cloneID == c).sum()), int((base.cloneID == c).sum()))
                   for c in shared if (sub.cloneID == c).sum() >= 10 and (base.cloneID == c).sum() >= 10]
            if per:
                diffs = np.array([p[1] for p in per]); w = np.array([min(p[2], p[3]) for p in per], float)
                wm = float((diffs * w).sum() / w.sum())
                bs = [float((diffs[i] * w[i]).sum() / w[i].sum())
                      for i in [rng_global.choice(len(diffs), len(diffs), True) for _ in range(2000)]]
                lo, hi = np.percentile(bs, [2.5, 97.5])
            else:
                wm = lo = hi = np.nan
            rows_h1.append(dict(mouse=m, site=s, n_site=len(sub), n_pt=len(base), raw_diff=raw,
                                within_clone_diff=wm, lo=float(lo), hi=float(hi), n_clones_paired=len(per),
                                mw_p=float(stats.mannwhitneyu(sub.e, base.e).pvalue),
                                frac_LH_site=float((sub.state4 == 2).mean()),
                                frac_LH_pt=float((base.state4 == 2).mean())))

        # ---- H3: clone PT mean vs clone metastatic mean
        for s in MET_SITES:
            sub = d[d.site == s]
            shared = [c for c in set(sub.cloneID) & set(base.cloneID)
                      if (sub.cloneID == c).sum() >= 20 and (base.cloneID == c).sum() >= 20]
            if len(shared) >= 4:
                x = np.array([base[base.cloneID == c].e.mean() for c in shared])
                y = np.array([sub[sub.cloneID == c].e.mean() for c in shared])
                r = stats.pearsonr(x, y)
                perm = [stats.pearsonr(x, rng_global.permutation(y)).statistic for _ in range(5000)]
                rows_h3.append(dict(mouse=m, site=s, n_clones=len(shared), r=float(r.statistic),
                                    p_param=float(r.pvalue),
                                    p_perm=float((np.abs(perm) >= abs(r.statistic)).mean()),
                                    clones=";".join(map(str, sorted(shared)))))

        # ---- H2/H4/H5: covariance decomposition per site
        for s in ["PT"] + MET_SITES + (["PTab", "Blood"] if m == 1 else []):
            sub = d[d.site == s]
            if len(sub) < 150:
                continue
            paths = list(sub.path.values)
            res = decompose(paths, sub.state2.values, sub.cloneID.values, n_boot=400, rng=rng_global)
            if res:
                rows_dec.append(dict(mouse=m, site=s, confirmatory=s in (["PT"] + MET_SITES), **res))
            cf = curve_frame(paths, sub.state2.values, sub.cloneID.values, rng=rng_global)
            if cf is not None:
                rows_curve += [dict(mouse=m, site=s, depth=int(a), clone_centred=float(b), site_centred=float(c),
                                    n_pairs=int(n)) for a, b, c, n in cf]
            # H2: continuous-score version, is |e_i - e_j| smaller for deeper MRCA
            idx = sub.index.values
            sm = sub.loc[idx]
            pick = rng_global.choice(len(sm), min(len(sm), 1200), replace=False)
            ee = sm.e.values[pick]; pp = [paths[i] for i in pick]; cc = sm.cloneID.values[pick]
            dd, ab = [], []
            for i in range(len(pick)):
                for j in range(i + 1, len(pick)):
                    if cc[i] == cc[j]:
                        dd.append(shared_depth(pp[i], pp[j])); ab.append(abs(ee[i] - ee[j]))
            if len(dd) > 200:
                dd = np.array(dd); ab = np.array(ab)
                r = stats.spearmanr(dd, ab)
                rows_h2.append(dict(mouse=m, site=s, n_pairs=len(dd), spearman=float(r.statistic),
                                    p=float(r.pvalue),
                                    mean_abs_shallow=float(ab[dd <= np.median(dd)].mean()),
                                    mean_abs_deep=float(ab[dd > np.median(dd)].mean())))

        # ---- exploratory: held-out programmes vs clone propensity
        for prog in HELDOUT:
            cl = d.groupby("cloneID").agg(n=("e", "size"), pr=(prog, "mean"), em=("e", "mean"))
            cl = cl[cl.n >= 30]
            if len(cl) >= 5:
                r = stats.spearmanr(cl.em, cl.pr)
                rows_prog.append(dict(mouse=m, prog=prog, n_clones=len(cl), rho=float(r.statistic),
                                      p=float(r.pvalue)))
        d.drop(columns=["path"]).to_csv(os.path.join(tab, f"m2_real_M{m}_cells.csv.gz"), index=False)

    pd.DataFrame(rows_h1).to_csv(os.path.join(tab, "m2_real_h1_enrichment.csv"), index=False)
    pd.DataFrame(rows_h3).to_csv(os.path.join(tab, "m2_real_h3_clone_corr.csv"), index=False)
    pd.DataFrame(rows_dec).to_csv(os.path.join(tab, "m2_real_decomposition.csv"), index=False)
    pd.DataFrame(rows_curve).to_csv(os.path.join(tab, "m2_real_curves.csv"), index=False)
    pd.DataFrame(rows_h2).to_csv(os.path.join(tab, "m2_real_h2_heritability.csv"), index=False)
    pd.DataFrame(rows_prog).to_csv(os.path.join(tab, "m2_real_programmes.csv"), index=False)
    with open(os.path.join(met, "m2_real_metrics.json"), "w") as fh:
        json.dump(M, fh, indent=1, default=float)
    print(json.dumps(M, indent=1, default=float))


if __name__ == "__main__":
    main(*sys.argv[1:5])
