"""E16: verify the capacity law against the ACTUAL geometry of the published lineage trees.

The law is derived for a balanced resolved tree. Real trees are unbalanced, so the covered band must be
recomputed from the empirical distribution of within-clone pair distances rather than from D_eff and D.

For each real clone:
    L_max_emp = the 90th percentile of pair lineage distance (robust to a few very deep pairs)
    L_min_emp = the 10th percentile, floored at one resolved branch
    covered k in [0.10 / L_max_emp, 1.0 / L_min_emp]
and this is compared with the idealised prediction [0.10 / D, D_eff / D].

Usage: python -I capacity_verification.py <tables_dir>
"""
import os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LO_KL, HI_KL = 0.10, 1.00


def main(tab):
    f = os.path.join(tab, "m3_depth_spectra.csv")
    fs = os.path.join(tab, "m3_depth_spectra_summary.csv")
    if not (os.path.exists(f) and os.path.exists(fs)):
        print("depth spectra missing"); return
    spec = pd.read_csv(f)          # system, dataset, clone, depth, frac
    summ = pd.read_csv(fs)
    rows = []
    for (sysname, ds, clone), g in spec.groupby(["system", "dataset", "clone"]):
        g = g.sort_values("depth")
        dmax = g.depth.max()
        # lineage distance from MRCA to each leaf, in resolved branches
        L = (dmax - g.depth.values).astype(float)
        w = g.frac.values / g.frac.values.sum()
        # empirical quantiles of the pair-distance distribution
        order = np.argsort(L); Ls = L[order]; ws = np.cumsum(w[order])
        q10 = float(np.interp(0.10, ws, Ls)); q90 = float(np.interp(0.90, ws, Ls))
        L_min = max(q10, 0.5); L_max = max(q90, L_min + 0.5)
        rows.append(dict(system=sysname, dataset=ds, clone=clone, depth_max=int(dmax),
                         L_min_emp=L_min, L_max_emp=L_max,
                         k_lo_emp=LO_KL / L_max, k_hi_emp=HI_KL / L_min,
                         decades_emp=np.log10((HI_KL / L_min) / (LO_KL / L_max))))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(tab, "e16_empirical_capacity.csv"), index=False)
    agg = d.groupby(["system", "dataset"]).agg(
        clones=("clone", "size"), depth_max=("depth_max", "median"),
        L_min=("L_min_emp", "median"), L_max=("L_max_emp", "median"),
        k_lo=("k_lo_emp", "median"), k_hi=("k_hi_emp", "median"),
        decades=("decades_emp", "median")).reset_index()
    # idealised prediction: D_eff = median max resolved depth, D unknown, so compare the DECADES only,
    # which the law says is 1 + log10(D_eff) and independent of D
    agg["decades_predicted"] = 1 + np.log10(agg.depth_max.clip(lower=1))
    agg["abs_error"] = (agg.decades - agg.decades_predicted).abs()
    agg.to_csv(os.path.join(tab, "e16_capacity_vs_prediction.csv"), index=False)
    print(agg.round(3).to_string(index=False))
    print("\nmedian |empirical - predicted| decades: %.3f" % agg.abs_error.median())
    print("\nfraction of clones whose covered band includes k = 1 per resolved branch: %.3f"
          % float((d.k_hi_emp >= 1.0).mean()))
    print("fraction covering k >= 2: %.3f" % float((d.k_hi_emp >= 2.0).mean()))


if __name__ == "__main__":
    main(sys.argv[1])
