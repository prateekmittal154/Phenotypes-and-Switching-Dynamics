"""Animal-level comparison of method conclusions on the Simeonov and Quinn datasets.

There is no ground truth here, so the question is not accuracy but agreement: do the
methods reach the same scientific conclusion from the same cells and the same trees, and
where they differ, which assumption is responsible?

Animals are kept separate throughout. Cells are not pooled across animals.

Usage: python -I analyse_real_data.py <benchdir>
"""
import os, sys, glob
import numpy as np
import pandas as pd

B = sys.argv[1]
MET = os.path.join(B, "results", "metrics")
TAB = os.path.join(B, "results", "tables")
os.makedirs(TAB, exist_ok=True)

comp = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(MET, "real_*_comp.csv")))
                  if "depth_sweep" not in f], ignore_index=True)
sd = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(MET, "real_*_sd.csv")))],
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

# SCOUT on the real trees was re-run separately after the collapse.singles fix
# (reports/deviations.md, B12). Its columns replace the failed ones from the first pass.
_sc = os.path.join(MET, "real_scout_all.csv")
if os.path.exists(_sc):
    sc = _norm(pd.read_csv(_sc))
    keep = ["rep_dir", "tree_variant", "clone", "SCOUT_alpha", "SCOUT_theta",
            "SCOUT_sigma2", "SCOUT_model", "SCOUT_note"]
    sc = sc[[c for c in keep if c in sc.columns]]
    df = df.drop(columns=[c for c in ("SCOUT_alpha", "SCOUT_theta", "SCOUT_sigma2",
                                      "SCOUT_model", "SCOUT_note") if c in df.columns])
    df = df.merge(sc, on=key, how="left")
    print("merged re-run SCOUT columns for",
          int(pd.to_numeric(df.SCOUT_alpha, errors="coerce").notna().sum()), "clones")

# rep_dir basename is "<dataset>__<clone>"; recover the animal from it
df["animal"] = df.rep_dir.apply(lambda p: os.path.basename(p).split("__")[0])
df["n_cells"] = pd.to_numeric(df["n_cells"], errors="coerce")
df.to_csv(os.path.join(MET, "real_clone_level_all.csv"), index=False)
print("real clone-level rows:", len(df))
print(df.groupby("animal").size().to_string())

ANIMAL_LABEL = {"simeonov_M1": "Simeonov M1", "simeonov_M2": "Simeonov M2",
                "quinn_5k": "Quinn 5k", "quinn_10k": "Quinn 10k",
                "quinn_30k": "Quinn 30k", "quinn_100k": "Quinn 100k"}

rows = []
for animal, g in df.groupby("animal"):
    n_clone = len(g)
    for meth, pcol, kcol in [("Study D", "StudyD_pi", "StudyD_k"),
                             ("PATH", "PATH_pi", "PATH_k"),
                             ("PASTRI", "PASTRI_pi", "PASTRI_k"),
                             ("SCOUT", "SCOUT_theta", "SCOUT_alpha")]:
        p = pd.to_numeric(g[pcol], errors="coerce")
        k = pd.to_numeric(g[kcol], errors="coerce")
        pf = p[np.isfinite(p)]
        kf = k[np.isfinite(k) & (k > 0)]
        # occupancy spread, scaled so the four methods are on a comparable 0-1 footing:
        # SCOUT's theta is on the continuous trait scale, the others on a probability scale.
        if len(pf) > 1:
            rng = pf.max() - pf.min()
            v_scaled = float(np.var(pf / rng)) if rng > 0 else 0.0
        else:
            v_scaled = np.nan
        rows.append(dict(
            animal=ANIMAL_LABEL.get(animal, animal), method=meth,
            n_clones=n_clone, n_with_occupancy=int(len(pf)), n_with_rate=int(len(kf)),
            occupancy_var_raw=float(np.var(pf)) if len(pf) > 1 else np.nan,
            V_pi_scaled=v_scaled,
            occupancy_range=float(pf.max() - pf.min()) if len(pf) > 1 else np.nan,
            rate_median=float(np.median(kf)) if len(kf) else np.nan,
            rate_iqr_log=float(np.subtract(*np.percentile(np.log(kf), [75, 25])))
            if len(kf) > 3 else np.nan,
            rate_var_log=float(np.var(np.log(kf))) if len(kf) > 1 else np.nan))
out = pd.DataFrame(rows)
out.to_csv(os.path.join(TAB, "TableB5_real_data.csv"), index=False)
print("\n=== Table B5: animal-level conclusions ===")
print(out.round(3).to_string(index=False))

# ---- cross-method agreement on which clones are high-occupancy, within each animal
agree = []
for animal, g in df.groupby("animal"):
    if len(g) < 4:
        continue
    cols = {"Study D": "StudyD_pi", "PATH": "PATH_pi", "PASTRI": "PASTRI_pi",
            "SCOUT": "SCOUT_theta"}
    ms = list(cols)
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            a = pd.to_numeric(g[cols[ms[i]]], errors="coerce")
            b = pd.to_numeric(g[cols[ms[j]]], errors="coerce")
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() < 4:
                continue
            r = float(pd.Series(a[ok]).corr(pd.Series(b[ok]), method="spearman"))
            agree.append(dict(animal=ANIMAL_LABEL.get(animal, animal),
                              pair=f"{ms[i]} vs {ms[j]}", quantity="occupancy",
                              n_clones=int(ok.sum()), spearman=r))
    cols = {"Study D": "StudyD_k", "PATH": "PATH_k", "PASTRI": "PASTRI_k",
            "SCOUT": "SCOUT_alpha"}
    ms = list(cols)
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            a = pd.to_numeric(g[cols[ms[i]]], errors="coerce")
            b = pd.to_numeric(g[cols[ms[j]]], errors="coerce")
            ok = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
            if ok.sum() < 4:
                continue
            r = float(pd.Series(np.log(a[ok])).corr(pd.Series(np.log(b[ok])), method="spearman"))
            agree.append(dict(animal=ANIMAL_LABEL.get(animal, animal),
                              pair=f"{ms[i]} vs {ms[j]}", quantity="rate",
                              n_clones=int(ok.sum()), spearman=r))
ag = pd.DataFrame(agree)
ag.to_csv(os.path.join(TAB, "TableB5b_cross_method_agreement.csv"), index=False)
print("\n=== Cross-method rank agreement, within animal ===")
if len(ag):
    piv = ag.pivot_table(index="pair", columns="quantity", values="spearman", aggfunc="mean")
    print(piv.round(3).to_string())
    print("\nmedian occupancy agreement:",
          round(float(ag[ag.quantity == 'occupancy'].spearman.median()), 3))
    print("median rate agreement:     ",
          round(float(ag[ag.quantity == 'rate'].spearman.median()), 3))
print("\nDone.")
