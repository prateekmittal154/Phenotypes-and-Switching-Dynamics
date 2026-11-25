"""Controls for the Simeonov re-analysis (deviation D5, added after H1 showed a within-clone reversal).

C1  Tissue-technical control. The within-clone primary-to-metastasis shift is recomputed for every held-out
    programme and for QC covariates (genes detected, mitochondrial and ribosomal fraction). If the EMT shift is
    biological it should not be accompanied by a shift of every unrelated programme in the same direction.
C2  Clone-stratified permutation test for H2 (the pre-registered test; the Spearman p-value over non-independent
    pairs is not valid).
C3  Decomposition of the between-clone component: how much of the population primary-to-metastasis difference is
    explained by the change in clone composition alone.
C4  Sign test over clones for the within-clone shift.
"""
import json, os, sys
import numpy as np, pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.biological_validation.extract_panel import HELDOUT
from src.lineage_inference.propensity import shared_depth

MET_SITES = ["Met", "Liver", "Lung"]
QC = ["nFeature_RNA", "percent_mt", "percent_ribo"]
rng = np.random.default_rng(20261127)


def main(tab, met, proc):
    rows, rows_perm, rows_dec, rows_sign = [], [], [], []
    for m in (1, 2):
        d = pd.read_csv(os.path.join(tab, f"m2_real_M{m}_cells.csv.gz"))
        base = d[d.site == "PT"]
        feats = ["e"] + list(HELDOUT) + [q for q in QC if q in d.columns]
        for s in MET_SITES:
            sub = d[d.site == s]
            if len(sub) < 50:
                continue
            shared = [c for c in set(sub.cloneID) & set(base.cloneID)
                      if (sub.cloneID == c).sum() >= 10 and (base.cloneID == c).sum() >= 10]
            if not shared:
                continue
            for f in feats:
                diffs = np.array([sub[sub.cloneID == c][f].mean() - base[base.cloneID == c][f].mean()
                                  for c in shared])
                w = np.array([min((sub.cloneID == c).sum(), (base.cloneID == c).sum()) for c in shared], float)
                # standardise by the pooled primary SD so programmes are comparable
                sd = base[f].std() or 1.0
                wm = float((diffs * w).sum() / w.sum() / sd)
                rows.append(dict(mouse=m, site=s, feature=f, n_clones=len(shared),
                                 within_clone_shift_sd=wm,
                                 raw_shift_sd=float((sub[f].mean() - base[f].mean()) / sd),
                                 n_neg=int((diffs < 0).sum())))
            # C4 sign test on the EMT score
            diffs = np.array([sub[sub.cloneID == c].e.mean() - base[base.cloneID == c].e.mean() for c in shared])
            rows_sign.append(dict(mouse=m, site=s, n_clones=len(shared), n_negative=int((diffs < 0).sum()),
                                  p_sign=float(stats.binomtest(int((diffs < 0).sum()), len(diffs), 0.5).pvalue)))
            # C3 composition vs within-clone decomposition of the population difference
            pc = base.cloneID.value_counts(normalize=True)
            qc_ = sub.cloneID.value_counts(normalize=True)
            cl = sorted(set(pc.index) | set(qc_.index))
            mu_pt = {c: base[base.cloneID == c].e.mean() for c in cl}
            mu_mt = {c: sub[sub.cloneID == c].e.mean() for c in cl}
            comp = sum(qc_.get(c, 0) * mu_pt.get(c, np.nan) for c in cl if np.isfinite(mu_pt.get(c, np.nan)))
            comp0 = sum(pc.get(c, 0) * mu_pt.get(c, np.nan) for c in cl if np.isfinite(mu_pt.get(c, np.nan)))
            within = sum(qc_.get(c, 0) * (mu_mt.get(c, np.nan) - mu_pt.get(c, np.nan))
                         for c in cl if np.isfinite(mu_mt.get(c, np.nan)) and np.isfinite(mu_pt.get(c, np.nan)))
            rows_dec.append(dict(mouse=m, site=s, total=float(sub.e.mean() - base.e.mean()),
                                 composition=float(comp - comp0), within_clone=float(within)))
        # C2 clone-stratified permutation for H2
        for s in ["PT"] + MET_SITES:
            sub = d[d.site == s]
            if len(sub) < 150:
                continue
            el = pd.read_csv(os.path.join(proc, f"M{m}/output-files/pdac_mouse{m}_full_edgelist.txt"), sep="\t")
            el["V1"] = el.V1.str.strip(); el["V2"] = el.V2.str.strip()
            par = dict(zip(el.V2, el.V1))
            def path(leaf):
                p, n = [], leaf
                while n in par:
                    n = par[n]
                    if n != "Cancer_root":
                        p.append(n)
                return tuple(reversed(p))
            pick = rng.choice(len(sub), min(len(sub), 1000), replace=False)
            sm = sub.iloc[pick]
            pp = [path(l) for l in sm.leaf]; ee = sm.e.values; cc = sm.cloneID.values
            dd, ab, cl = [], [], []
            for i in range(len(sm)):
                for j in range(i + 1, len(sm)):
                    if cc[i] == cc[j]:
                        dd.append(shared_depth(pp[i], pp[j])); ab.append(abs(ee[i] - ee[j])); cl.append(cc[i])
            if len(dd) < 200:
                continue
            dd, ab, cl = np.array(dd), np.array(ab), np.array(cl)
            obs = stats.spearmanr(dd, ab).statistic
            null = []
            for _ in range(1000):
                ee2 = ee.copy()
                for c in np.unique(cc):                      # permute states within clone, keeping the tree
                    k = np.where(cc == c)[0]
                    ee2[k] = rng.permutation(ee2[k])
                ab2 = np.abs(ee2[[0]] * 0)  # placeholder
                ab2 = np.array([abs(ee2[i] - ee2[j]) for i in range(len(sm))
                                for j in range(i + 1, len(sm)) if cc[i] == cc[j]])
                null.append(stats.spearmanr(dd, ab2).statistic)
            null = np.array(null)
            rows_perm.append(dict(mouse=m, site=s, n_pairs=len(dd), obs_spearman=float(obs),
                                  null_mean=float(null.mean()), null_sd=float(null.std()),
                                  p_perm=float((null <= obs).mean())))
            print("perm done", m, s, flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(tab, "m2_real_c1_technical_control.csv"), index=False)
    pd.DataFrame(rows_perm).to_csv(os.path.join(tab, "m2_real_c2_h2_permutation.csv"), index=False)
    pd.DataFrame(rows_dec).to_csv(os.path.join(tab, "m2_real_c3_composition.csv"), index=False)
    pd.DataFrame(rows_sign).to_csv(os.path.join(tab, "m2_real_c4_signtest.csv"), index=False)
    print("controls done")


if __name__ == "__main__":
    main(*sys.argv[1:4])
