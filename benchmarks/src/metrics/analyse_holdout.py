"""Estimand group D: held-out prediction, the one axis all four methods share.

Each method's per-clone central tendency predicts the withheld cells of that clone.
Binary methods predict the state with their fitted occupancy; SCOUT predicts the
continuous trait with its fitted OU optimum. The two are made comparable by a skill
score against the same baseline, the grand mean of the training cells:

    skill = 1 - MSE(method) / MSE(grand mean)

A skill of 0 means a method adds nothing over ignoring clone identity.

Usage: python -I analyse_holdout.py <benchdir>
"""
import os, sys, csv, glob
import numpy as np
import pandas as pd

B = sys.argv[1]
MET = os.path.join(B, "results", "metrics")
TAB = os.path.join(B, "results", "tables")
os.makedirs(TAB, exist_ok=True)


def _norm(t):
    if "scenario" in t.columns:
        t = t.drop(columns=["scenario"])
    t["tree_variant"] = t["tree_variant"].astype(str).str.strip().str.lower()
    t["clone"] = t["clone"].astype(str)
    t["rep_dir"] = t["rep_dir"].astype(str).str.rstrip("/")
    return t


comp = _norm(pd.read_csv(os.path.join(MET, "holdout_comp.csv")))
sd = _norm(pd.read_csv(os.path.join(MET, "holdout_sd.csv")))
df = comp.merge(sd, on=["rep_dir", "tree_variant", "clone"], how="inner", suffixes=("", "_sd"))
print("held-out clone rows:", len(df))

rows = []
for rep, g in df.groupby("rep_dir"):
    ho = os.path.join(rep, "heldout.csv")
    tr = os.path.join(rep, "states.csv")
    if not (os.path.exists(ho) and os.path.exists(tr)):
        continue
    test = list(csv.DictReader(open(ho)))
    train = list(csv.DictReader(open(tr)))
    if not test:
        continue
    base_state = np.mean([int(r["state"]) for r in train])
    base_trait = np.mean([float(r["trait"]) for r in train])
    pi_map = {m: dict(zip(g.clone, pd.to_numeric(g[c], errors="coerce")))
              for m, c in [("Study D", "StudyD_pi"), ("PATH", "PATH_pi"),
                           ("PASTRI", "PASTRI_pi")]}
    th_map = dict(zip(g.clone, pd.to_numeric(g["SCOUT_theta"], errors="coerce")))

    y_state = np.array([int(r["state"]) for r in test], float)
    y_trait = np.array([float(r["trait"]) for r in test], float)
    clones = [r["clone"] for r in test]
    mse_base_s = float(np.mean((y_state - base_state) ** 2))
    mse_base_t = float(np.mean((y_trait - base_trait) ** 2))

    for m in ["Study D", "PATH", "PASTRI"]:
        pred = np.array([pi_map[m].get(c, np.nan) for c in clones], float)
        ok = np.isfinite(pred)
        if ok.sum() < 5:
            continue
        mse = float(np.mean((y_state[ok] - pred[ok]) ** 2))
        rows.append(dict(rep_dir=rep, method=m, target="binary state",
                         n_heldout=int(ok.sum()), mse=mse, mse_baseline=mse_base_s,
                         skill=1 - mse / mse_base_s if mse_base_s > 0 else np.nan,
                         mae=float(np.mean(np.abs(y_state[ok] - pred[ok])))))
    pred = np.array([th_map.get(c, np.nan) for c in clones], float)
    ok = np.isfinite(pred)
    if ok.sum() >= 5:
        mse = float(np.mean((y_trait[ok] - pred[ok]) ** 2))
        rows.append(dict(rep_dir=rep, method="SCOUT", target="continuous trait",
                         n_heldout=int(ok.sum()), mse=mse, mse_baseline=mse_base_t,
                         skill=1 - mse / mse_base_t if mse_base_t > 0 else np.nan,
                         mae=float(np.mean(np.abs(y_trait[ok] - pred[ok])))))

out = pd.DataFrame(rows)
out.to_csv(os.path.join(MET, "holdout_per_rep.csv"), index=False)


def boot(x, B_=2000, seed=0):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    if len(x) < 3:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    s = [np.median(rng.choice(x, len(x), replace=True)) for _ in range(B_)]
    return float(np.percentile(s, 2.5)), float(np.percentile(s, 97.5))


summ = []
for (m, t), g in out.groupby(["method", "target"]):
    lo, hi = boot(g.skill.values)
    summ.append(dict(method=m, target=t, n_reps=len(g),
                     median_skill=float(np.nanmedian(g.skill)), skill_lo=lo, skill_hi=hi,
                     median_mae=float(np.nanmedian(g.mae))))
s = pd.DataFrame(summ).sort_values("median_skill", ascending=False)
s.to_csv(os.path.join(TAB, "TableB3e_holdout_prediction.csv"), index=False)
print("\n=== Held-out prediction skill (1 = perfect, 0 = no better than the grand mean) ===")
print(s.round(4).to_string(index=False))
