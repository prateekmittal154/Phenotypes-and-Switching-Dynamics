"""Technical control for the Quinn within-clone shifts, mirroring the Simeonov control (Milestone 2, D5).

If the within-clone primary-to-metastasis shift in the cell-cycle state is biological, it should not be
accompanied by shifts of the same size in quality covariates.
Usage: python -I quinn_controls.py <states_dir> <tables_dir>
"""
import os, sys
import numpy as np, pandas as pd

PRIMARY = ("LL", "RL", "RE", "RW", "NE", "NW", "SE", "SW")
FEATS = ["cc_score", "pc1", "mt_frac", "n_genes", "lib"]


def main(sdir, tab):
    rows = []
    for mouse in ("5k", "30k", "10k", "100k"):
        f = os.path.join(sdir, f"quinn_cells_{mouse}.csv.gz")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        d["site"] = d.sampleID.astype(str).str.replace("^30k", "", regex=True)
        base = d[d.site.isin(PRIMARY)]
        if len(base) < 50:
            continue
        for s in sorted(set(d.site) - set(PRIMARY)):
            sub = d[d.site == s]
            if len(sub) < 40:
                continue
            shared = [c for c in set(sub.LineageGroup) & set(base.LineageGroup)
                      if (sub.LineageGroup == c).sum() >= 10 and (base.LineageGroup == c).sum() >= 10]
            if not shared:
                continue
            w = np.array([min((sub.LineageGroup == c).sum(), (base.LineageGroup == c).sum())
                          for c in shared], float)
            for feat in FEATS:
                if feat not in d.columns:
                    continue
                diffs = np.array([sub[sub.LineageGroup == c][feat].mean() -
                                  base[base.LineageGroup == c][feat].mean() for c in shared])
                sd = base[feat].std() or 1.0
                rows.append(dict(mouse=mouse, site=s, feature=feat, n_clones=len(shared),
                                 within_clone_shift_sd=float((diffs * w).sum() / w.sum() / sd),
                                 raw_shift_sd=float((sub[feat].mean() - base[feat].mean()) / sd)))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(tab, "m4_quinn_technical_control.csv"), index=False)
    print(out.pivot_table(index=["mouse", "site"], columns="feature",
                          values="within_clone_shift_sd").round(3).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:3])
