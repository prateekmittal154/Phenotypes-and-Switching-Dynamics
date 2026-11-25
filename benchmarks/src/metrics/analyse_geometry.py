"""Benchmark II: does the theoretical information window predict real method performance?

Across five tree geometries the true switching rate is identical (k = 0.3 for every clone,
so the true V_log k is exactly zero). Anything a method reports as rate heterogeneity is
therefore produced by geometry alone. The independent variable is the fraction of sampled
within-clone pairs whose kL falls in the informative band [0.10, 1.0].

Usage: python -I analyse_geometry.py <benchdir>
"""
import os, sys, json, glob
import numpy as np
import pandas as pd

B = sys.argv[1]
MET = os.path.join(B, "results", "metrics")
TAB = os.path.join(B, "results", "tables")
os.makedirs(TAB, exist_ok=True)

comp = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(MET, "geom_*_comp.csv")))
                  if "depth_sweep" not in f], ignore_index=True)
sd = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(MET, "geom_*_sd.csv")))],
               ignore_index=True)
def _norm(t):
    """`true` is read back as a bool from one writer and a string from the other."""
    if "scenario" in t.columns:
        t = t.drop(columns=["scenario"])
    t["tree_variant"] = (t["tree_variant"].astype(str).str.strip().str.lower()
                         .replace({"true": "true", "unit": "unit"}))
    t["clone"] = t["clone"].astype(str)
    t["rep_dir"] = t["rep_dir"].astype(str).str.rstrip("/")
    return t

# The two adapters label `scenario` differently, so join on the replicate directory.
comp, sd = _norm(comp), _norm(sd)
key = ["rep_dir", "tree_variant", "clone"]
df = comp.merge(sd, on=key, how="inner", suffixes=("", "_sd"))
df["geometry"] = df.rep_dir.apply(lambda p: "_".join(os.path.basename(p).split("_")[:2]))

geo = {}
for d in df.rep_dir.unique():
    tj = os.path.join(d, "truth.json")
    if os.path.exists(tj):
        geo[d] = json.load(open(tj))
df["frac_in_window"] = df.rep_dir.map(lambda d: geo.get(d, {}).get("frac_in_window", np.nan))
df["frac_close"] = df.rep_dir.map(lambda d: geo.get(d, {}).get("frac_close", np.nan))
df["mrca_median"] = df.rep_dir.map(lambda d: geo.get(d, {}).get("mrca_median", np.nan))
df["true_k"] = 0.3
df.to_csv(os.path.join(MET, "geometry_clone_level.csv"), index=False)
print("geometry clone rows:", len(df))

rows = []
for (g, rep), gg in df.groupby(["geometry", "rep_dir"]):
    base = dict(geometry=g, rep_dir=rep,
                frac_in_window=gg.frac_in_window.iloc[0],
                frac_close=gg.frac_close.iloc[0],
                mrca_median=gg.mrca_median.iloc[0])
    for meth, kcol in [("Study D", "StudyD_k"), ("PATH", "PATH_k"),
                       ("PASTRI", "PASTRI_k"), ("SCOUT", "SCOUT_alpha")]:
        k = pd.to_numeric(gg[kcol], errors="coerce")
        k = k[np.isfinite(k) & (k > 0)]
        rows.append(dict(**base, method=meth,
                         # true V_log k is 0 by construction: this is spurious heterogeneity
                         spurious_V_logk=float(np.var(np.log(k))) if len(k) > 1 else np.nan,
                         median_log_bias=float(np.median(np.log(k) - np.log(0.3)))
                         if len(k) else np.nan,
                         n_clones=len(k)))
per_rep = pd.DataFrame(rows)
per_rep.to_csv(os.path.join(MET, "geometry_per_rep.csv"), index=False)

summ = per_rep.groupby(["geometry", "method"]).agg(
    n_reps=("spurious_V_logk", "size"),
    frac_in_window=("frac_in_window", "mean"),
    frac_close=("frac_close", "mean"),
    mrca_median=("mrca_median", "mean"),
    spurious_V_logk=("spurious_V_logk", "median"),
    median_log_bias=("median_log_bias", "median")).reset_index()
summ.to_csv(os.path.join(TAB, "TableB3c_geometry.csv"), index=False)

print("\n=== Spurious rate heterogeneity by tree geometry (true V_log k = 0) ===")
print(summ.pivot_table(index=["geometry", "frac_in_window"], columns="method",
                       values="spurious_V_logk").round(3).to_string())
print("\n=== Median log bias in the rate (true k = 0.3 everywhere) ===")
print(summ.pivot_table(index=["geometry", "frac_in_window"], columns="method",
                       values="median_log_bias").round(3).to_string())

# does informativeness predict performance?
print("\n=== Does the information window predict performance? ===")
pred = []
for meth, g in per_rep.groupby("method"):
    ok = np.isfinite(g.frac_in_window) & np.isfinite(g.spurious_V_logk)
    if ok.sum() < 5:
        continue
    r_sp = float(pd.Series(g.frac_in_window[ok]).corr(pd.Series(g.spurious_V_logk[ok]),
                                                      method="spearman"))
    ok2 = np.isfinite(g.frac_in_window) & np.isfinite(g.median_log_bias)
    r_b = float(pd.Series(g.frac_in_window[ok2]).corr(pd.Series(np.abs(g.median_log_bias[ok2])),
                                                      method="spearman"))
    pred.append(dict(method=meth, n_reps=int(ok.sum()),
                     spearman_window_vs_spurious_Vlogk=r_sp,
                     spearman_window_vs_abs_bias=r_b))
pdf = pd.DataFrame(pred)
pdf.to_csv(os.path.join(TAB, "TableB3d_window_predicts.csv"), index=False)
print(pdf.round(3).to_string(index=False))
print("\nNegative correlations mean more informative trees give less spurious "
      "heterogeneity and smaller bias, i.e. the window predicts performance.")
