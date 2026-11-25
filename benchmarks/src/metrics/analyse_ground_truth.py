"""Combine per-clone estimates from all four methods and score them against the truth.

Scoring follows protocols/estimand_matrix.md. Methods are compared only inside the
admissible groups defined there; no single leaderboard is produced.

Usage: python -I analyse_ground_truth.py <benchdir>
"""
import os, sys, json, glob
import numpy as np
import pandas as pd

B = sys.argv[1]
MET = os.path.join(B, "results", "metrics")
TAB = os.path.join(B, "results", "tables")
os.makedirs(TAB, exist_ok=True)

SCEN_LABEL = {
    "S1": "stable-state selection",
    "S2": "equal occupancy, rate varies",
    "S3": "occupancy varies, equal rate",
    "S4": "joint heterogeneity",
    "S5": "environmental induction",
    "S6": "founder bottleneck",
    "S7": "unequal clonal expansion",
    "S8": "observation misspecification",
}


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def boot_ci(x, fn=np.median, B_=2000, seed=0):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    if len(x) < 3:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    s = [fn(rng.choice(x, len(x), replace=True)) for _ in range(B_)]
    return (float(np.percentile(s, 2.5)), float(np.percentile(s, 97.5)))


# ---------------------------------------------------------------- load
comp = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(MET, "gt_S*_*.csv")))
                  if "depth_sweep" not in f], ignore_index=True)
sd = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(MET, "sd_S*_*.csv")))],
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
df["scenario"] = df.rep_dir.apply(lambda p: os.path.basename(p).split("_")[0])
print("merged clone-level rows:", len(df),
      "| scenarios:", sorted(df.scenario.unique()),
      "| variants:", sorted(df.tree_variant.unique()))

# ---------------------------------------------------------------- attach truth
truth_rows = []
for d in sorted(df["rep_dir"].unique()):
    tj = os.path.join(d, "truth.json")
    if not os.path.exists(tj):
        continue
    t = json.load(open(tj))
    for c in t["clones"]:
        truth_rows.append(dict(rep_dir=d, clone=f"C{c['clone']}", true_pi=c["pi"],
                               true_k=c["k"], true_theta=c["theta"], true_alpha=c["alpha"],
                               true_V_pi=t["V_pi"], true_V_logk=t["V_logk"]))
tr = pd.DataFrame(truth_rows)
df = df.merge(tr, on=["rep_dir", "clone"], how="left")
df.to_csv(os.path.join(MET, "clone_level_all.csv"), index=False)
print("wrote clone_level_all.csv")

# ---------------------------------------------------------------- Table B3/B4 helpers
def summarise_clone_level(d, variant):
    """Per-scenario accuracy of each method's clone-level estimates."""
    out = []
    for scen, g in d[d.tree_variant == variant].groupby("scenario"):
        n_clone = len(g)
        # --- occupancy (Group B): correlation with truth and absolute error
        for meth, col in [("Study D", "StudyD_pi"), ("PATH", "PATH_pi"), ("PASTRI", "PASTRI_pi")]:
            m = g[col].astype(float)
            ok = np.isfinite(m) & np.isfinite(g.true_pi)
            r = np.corrcoef(m[ok], g.true_pi[ok])[0, 1] if ok.sum() > 3 else np.nan
            mae = float(np.mean(np.abs(m[ok] - g.true_pi[ok]))) if ok.sum() else np.nan
            lo, hi = boot_ci(np.abs(m[ok] - g.true_pi[ok]), np.mean)
            out.append(dict(scenario=scen, label=SCEN_LABEL[scen], quantity="occupancy",
                            method=meth, n=int(ok.sum()), n_possible=n_clone,
                            returned=float(ok.sum() / n_clone),
                            corr_with_truth=r, mean_abs_err=mae,
                            mae_lo=lo, mae_hi=hi))
        # SCOUT's occupancy analogue is the OU optimum theta, on the continuous trait
        m = g["SCOUT_theta"].astype(float)
        ok = np.isfinite(m) & np.isfinite(g.true_theta)
        r = np.corrcoef(m[ok], g.true_theta[ok])[0, 1] if ok.sum() > 3 else np.nan
        out.append(dict(scenario=scen, label=SCEN_LABEL[scen], quantity="occupancy",
                        method="SCOUT (theta)", n=int(ok.sum()), n_possible=n_clone,
                        returned=float(ok.sum() / n_clone), corr_with_truth=r,
                        mean_abs_err=np.nan, mae_lo=np.nan, mae_hi=np.nan))

        # --- switching rate (Group A): log-scale bias against the true k
        for meth, col in [("Study D", "StudyD_k"), ("PATH", "PATH_k"), ("PASTRI", "PASTRI_k")]:
            m = g[col].astype(float)
            ok = np.isfinite(m) & (m > 0) & np.isfinite(g.true_k)
            lr = np.log(m[ok]) - np.log(g.true_k[ok])
            med = float(np.median(lr)) if len(lr) else np.nan
            lo, hi = boot_ci(lr, np.median)
            out.append(dict(scenario=scen, label=SCEN_LABEL[scen], quantity="switching rate",
                            method=meth, n=int(ok.sum()), n_possible=n_clone,
                            returned=float(ok.sum() / n_clone),
                            corr_with_truth=(np.corrcoef(np.log(m[ok]), np.log(g.true_k[ok]))[0, 1]
                                             if ok.sum() > 3 and g.true_k[ok].std() > 0 else np.nan),
                            median_log_bias=med, bias_lo=lo, bias_hi=hi))
        m = g["SCOUT_alpha"].astype(float)
        ok = np.isfinite(m) & (m > 0) & np.isfinite(g.true_alpha)
        lr = np.log(m[ok]) - np.log(g.true_alpha[ok])
        out.append(dict(scenario=scen, label=SCEN_LABEL[scen], quantity="decay rate",
                        method="SCOUT (alpha)", n=int(ok.sum()), n_possible=n_clone,
                        returned=float(ok.sum() / n_clone),
                        median_log_bias=float(np.median(lr)) if len(lr) else np.nan,
                        corr_with_truth=(np.corrcoef(np.log(m[ok]), np.log(g.true_alpha[ok]))[0, 1]
                                         if ok.sum() > 3 and g.true_alpha[ok].std() > 0 else np.nan)))
    return pd.DataFrame(out)


def variance_components(d, variant):
    """Per-replicate across-clone variance components, the quantity the paper reports."""
    rows = []
    for (scen, rep), g in d[d.tree_variant == variant].groupby(["scenario", "rep_dir"]):
        base = dict(scenario=scen, label=SCEN_LABEL[scen], rep_dir=rep,
                    true_V_pi=g.true_V_pi.iloc[0], true_V_logk=g.true_V_logk.iloc[0])
        for meth, pcol, kcol in [("Study D", "StudyD_pi", "StudyD_k"),
                                 ("PATH", "PATH_pi", "PATH_k"),
                                 ("PASTRI", "PASTRI_pi", "PASTRI_k")]:
            p = pd.to_numeric(g[pcol], errors="coerce")
            k = pd.to_numeric(g[kcol], errors="coerce")
            p = p[np.isfinite(p)]
            k = k[np.isfinite(k) & (k > 0)]
            rows.append(dict(**base, method=meth,
                             V_pi=float(np.var(p)) if len(p) > 1 else np.nan,
                             V_logk=float(np.var(np.log(k))) if len(k) > 1 else np.nan,
                             n_clones_pi=len(p), n_clones_k=len(k)))
        th = pd.to_numeric(g["SCOUT_theta"], errors="coerce"); th = th[np.isfinite(th)]
        al = pd.to_numeric(g["SCOUT_alpha"], errors="coerce"); al = al[np.isfinite(al) & (al > 0)]
        rows.append(dict(**base, method="SCOUT",
                         V_pi=float(np.var(th)) if len(th) > 1 else np.nan,
                         V_logk=float(np.var(np.log(al))) if len(al) > 1 else np.nan,
                         n_clones_pi=len(th), n_clones_k=len(al)))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- Table B3 (true tree)
b3 = summarise_clone_level(df, "true")
b3.to_csv(os.path.join(TAB, "TableB3_clone_level_true_tree.csv"), index=False)
print("\n=== Table B3: clone-level accuracy, true branch lengths ===")
print(b3[b3.quantity == "occupancy"][["scenario", "method", "returned", "corr_with_truth",
                                      "mean_abs_err"]].round(3).to_string(index=False))
print()
print(b3[b3.quantity.isin(["switching rate", "decay rate"])][
    ["scenario", "method", "returned", "median_log_bias"]].round(3).to_string(index=False))

# ---------------------------------------------------------------- variance components
vc_true = variance_components(df, "true")
vc_unit = variance_components(df, "unit")
vc_true["tree_variant"] = "true"
vc_unit["tree_variant"] = "unit"
vc = pd.concat([vc_true, vc_unit], ignore_index=True)
vc.to_csv(os.path.join(MET, "variance_components_by_rep.csv"), index=False)

# does each method detect occupancy heterogeneity where it exists, and not where it doesn't?
summ = []
for (scen, meth, variant), g in vc.groupby(["scenario", "method", "tree_variant"]):
    summ.append(dict(scenario=scen, label=SCEN_LABEL[scen], method=meth, tree_variant=variant,
                     n_reps=len(g),
                     true_V_pi=float(np.nanmedian(g.true_V_pi)),
                     est_V_pi=float(np.nanmedian(g.V_pi)),
                     true_V_logk=float(np.nanmedian(g.true_V_logk)),
                     est_V_logk=float(np.nanmedian(g.V_logk))))
vs = pd.DataFrame(summ)
vs.to_csv(os.path.join(TAB, "TableB3b_variance_components.csv"), index=False)
print("\n=== Variance components (median over replicates), true branch lengths ===")
print(vs[vs.tree_variant == "true"][["scenario", "method", "true_V_pi", "est_V_pi",
                                     "true_V_logk", "est_V_logk"]].round(4).to_string(index=False))

# ---------------------------------------------------------------- Table B4: clock sensitivity
rows = []
piv = df.pivot_table(index=["rep_dir", "clone", "scenario"], columns="tree_variant",
                     values=["StudyD_k", "PATH_k", "PASTRI_k", "SCOUT_alpha",
                             "StudyD_pi", "PATH_pi", "PASTRI_pi"], aggfunc="first")
for meth, col in [("Study D", "StudyD_k"), ("PATH", "PATH_k"),
                  ("PASTRI", "PASTRI_k"), ("SCOUT", "SCOUT_alpha")]:
    try:
        a = pd.to_numeric(piv[(col, "true")], errors="coerce")
        b = pd.to_numeric(piv[(col, "unit")], errors="coerce")
    except KeyError:
        continue
    ok = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    lr = np.log(b[ok]) - np.log(a[ok])
    lo, hi = boot_ci(lr, np.median)
    rows.append(dict(method=meth, quantity="rate / decay", n=int(ok.sum()),
                     median_log_ratio_unit_vs_true=float(np.median(lr)) if len(lr) else np.nan,
                     ci_lo=lo, ci_hi=hi,
                     frac_identical=float(np.mean(np.isclose(a[ok], b[ok], rtol=1e-6)))))
for meth, col in [("Study D", "StudyD_pi"), ("PATH", "PATH_pi"), ("PASTRI", "PASTRI_pi")]:
    try:
        a = pd.to_numeric(piv[(col, "true")], errors="coerce")
        b = pd.to_numeric(piv[(col, "unit")], errors="coerce")
    except KeyError:
        continue
    ok = np.isfinite(a) & np.isfinite(b)
    rows.append(dict(method=meth, quantity="occupancy", n=int(ok.sum()),
                     median_abs_diff_unit_vs_true=float(np.median(np.abs(b[ok] - a[ok]))),
                     frac_identical=float(np.mean(np.isclose(a[ok], b[ok], rtol=1e-6)))))
b4 = pd.DataFrame(rows)
b4.to_csv(os.path.join(TAB, "TableB4_clock_sensitivity.csv"), index=False)
print("\n=== Table B4: effect of replacing true branch lengths with unit edges ===")
print(b4.round(4).to_string(index=False))

# ---------------------------------------------------------------- failure / boundary rates
fr = []
for (scen, variant), g in df.groupby(["scenario", "tree_variant"]):
    n = len(g)
    for meth, col in [("Study D", "StudyD_k"), ("PATH", "PATH_k"), ("PASTRI", "PASTRI_k"),
                      ("SCOUT", "SCOUT_alpha")]:
        v = pd.to_numeric(g[col], errors="coerce")
        nonfinite = int((~np.isfinite(v)).sum())
        lo, hi = wilson(nonfinite, n)
        fr.append(dict(scenario=scen, label=SCEN_LABEL[scen], tree_variant=variant, method=meth,
                       n=n, n_no_estimate=nonfinite, frac_no_estimate=nonfinite / n,
                       ci_lo=lo, ci_hi=hi))
    sb = pd.to_numeric(g.get("StudyD_k_at_bound", pd.Series(dtype=float)), errors="coerce")
    fr.append(dict(scenario=scen, label=SCEN_LABEL[scen], tree_variant=variant,
                   method="Study D (k at bound)", n=n,
                   n_no_estimate=int(np.nansum(sb)),
                   frac_no_estimate=float(np.nansum(sb)) / n, ci_lo=np.nan, ci_hi=np.nan))
fdf = pd.DataFrame(fr)
fdf.to_csv(os.path.join(TAB, "TableB4b_no_estimate_rates.csv"), index=False)
print("\n=== Fraction of clones where no finite estimate was returned (true tree) ===")
print(fdf[fdf.tree_variant == "true"].pivot_table(index="scenario", columns="method",
      values="frac_no_estimate").round(3).to_string())
print("\nDone.")
