"""Benchmark IV: what tree reconstruction costs each method.

Five trees per replicate over identical cells and states. Any difference in a method's
output is therefore attributable to the tree alone. Truth is k = 0.3 for every clone.

Usage: python -I analyse_reconstruction.py <benchdir>
"""
import os, sys, json, glob
import numpy as np
import pandas as pd

B = sys.argv[1]
MET = os.path.join(B, "results", "metrics")
TAB = os.path.join(B, "results", "tables")
os.makedirs(TAB, exist_ok=True)
VARIANTS = ["true", "unit", "greedy", "nj", "upgma"]


def _norm(t):
    if "scenario" in t.columns:
        t = t.drop(columns=["scenario"])
    t["tree_variant"] = t["tree_variant"].astype(str).str.strip().str.lower()
    t["clone"] = t["clone"].astype(str)
    t["rep_dir"] = t["rep_dir"].astype(str).str.rstrip("/")
    return t


comp = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(MET, "recon_*_comp.csv")))
                  if "depth_sweep" not in f], ignore_index=True)
sd = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(MET, "recon_*_sd.csv")))],
               ignore_index=True)
comp, sd = _norm(comp), _norm(sd)
df = comp.merge(sd, on=["rep_dir", "tree_variant", "clone"], how="inner", suffixes=("", "_sd"))
print("reconstruction clone rows:", len(df), "| variants:", sorted(df.tree_variant.unique()))

sat = {}
for d in df.rep_dir.unique():
    tj = os.path.join(d, "truth.json")
    if os.path.exists(tj):
        sat[d] = json.load(open(tj)).get("saturation", np.nan)
df["saturation"] = df.rep_dir.map(sat)
truth_pi = {}
for d in df.rep_dir.unique():
    tj = os.path.join(d, "truth.json")
    if os.path.exists(tj):
        for c in json.load(open(tj))["clones"]:
            truth_pi[(d, f"C{c['clone']}")] = c["pi"]
df["true_pi"] = [truth_pi.get((r, c), np.nan) for r, c in zip(df.rep_dir, df.clone)]
df["true_k"] = 0.3
df.to_csv(os.path.join(MET, "reconstruction_clone_level.csv"), index=False)

METHODS = [("Study D", "StudyD_k", "StudyD_pi"), ("PATH", "PATH_k", "PATH_pi"),
           ("PASTRI", "PASTRI_k", "PASTRI_pi"), ("SCOUT", "SCOUT_alpha", "SCOUT_theta")]

rows = []
for v in VARIANTS:
    g = df[df.tree_variant == v]
    if not len(g):
        continue
    for meth, kcol, pcol in METHODS:
        k = pd.to_numeric(g[kcol], errors="coerce")
        kf = k[np.isfinite(k) & (k > 0)]
        p = pd.to_numeric(g[pcol], errors="coerce")
        ok = np.isfinite(p) & np.isfinite(g.true_pi)
        r = float(pd.Series(p[ok]).corr(pd.Series(g.true_pi[ok]))) if ok.sum() > 3 else np.nan
        # spurious rate heterogeneity, per replicate then median
        sp = []
        for rep, gg in g.groupby("rep_dir"):
            kk = pd.to_numeric(gg[kcol], errors="coerce")
            kk = kk[np.isfinite(kk) & (kk > 0)]
            if len(kk) > 1:
                sp.append(float(np.var(np.log(kk))))
        rows.append(dict(tree_variant=v, method=meth, n=len(g),
                         frac_no_estimate=float((~np.isfinite(k)).mean()),
                         median_log_bias=float(np.median(np.log(kf) - np.log(0.3)))
                         if len(kf) else np.nan,
                         occupancy_corr=r,
                         spurious_V_logk=float(np.median(sp)) if sp else np.nan,
                         median_secs=float(pd.to_numeric(
                             g.get(meth.split()[0].upper() + "_secs",
                                   g.get("StudyD_secs")), errors="coerce").median())
                         if True else np.nan))
out = pd.DataFrame(rows)
out.to_csv(os.path.join(TAB, "TableB4c_reconstruction.csv"), index=False)

print("\n=== Occupancy correlation with truth, by tree ===")
print(out.pivot_table(index="tree_variant", columns="method",
                      values="occupancy_corr").reindex(VARIANTS).round(3).to_string())
print("\n=== Median log bias in the rate (true k = 0.3), by tree ===")
print(out.pivot_table(index="tree_variant", columns="method",
                      values="median_log_bias").reindex(VARIANTS).round(3).to_string())
print("\n=== Spurious rate heterogeneity (true V_log k = 0), by tree ===")
print(out.pivot_table(index="tree_variant", columns="method",
                      values="spurious_V_logk").reindex(VARIANTS).round(3).to_string())

# stability: how far each method's own estimate moves from the true-tree answer
stab = []
piv = df.pivot_table(index=["rep_dir", "clone"], columns="tree_variant",
                     values=["StudyD_k", "PATH_k", "PASTRI_k", "SCOUT_alpha",
                             "StudyD_pi", "PATH_pi", "PASTRI_pi", "SCOUT_theta"],
                     aggfunc="first")
for meth, kcol, pcol in METHODS:
    for v in VARIANTS:
        if v == "true":
            continue
        try:
            a = pd.to_numeric(piv[(kcol, "true")], errors="coerce")
            b = pd.to_numeric(piv[(kcol, v)], errors="coerce")
        except KeyError:
            continue
        ok = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
        try:
            ap = pd.to_numeric(piv[(pcol, "true")], errors="coerce")
            bp = pd.to_numeric(piv[(pcol, v)], errors="coerce")
            okp = np.isfinite(ap) & np.isfinite(bp)
            dp = float(np.median(np.abs(bp[okp] - ap[okp]))) if okp.sum() else np.nan
        except KeyError:
            dp = np.nan
        stab.append(dict(method=meth, tree_variant=v, n=int(ok.sum()),
                         median_abs_log_shift_rate=float(np.median(np.abs(np.log(b[ok]) -
                                                                          np.log(a[ok]))))
                         if ok.sum() else np.nan,
                         median_abs_shift_occupancy=dp))
st = pd.DataFrame(stab)
st.to_csv(os.path.join(TAB, "TableB4d_reconstruction_stability.csv"), index=False)
print("\n=== Movement away from the true-tree answer ===")
print(st.pivot_table(index="tree_variant", columns="method",
                     values="median_abs_log_shift_rate").round(3).to_string())
print("\n(occupancy)")
print(st.pivot_table(index="tree_variant", columns="method",
                     values="median_abs_shift_occupancy").round(4).to_string())
print("\nmean recorder saturation:", round(float(np.nanmean(list(sat.values()))), 3))
