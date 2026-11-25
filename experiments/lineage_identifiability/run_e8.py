"""E8: can the MRCA-depth covariance plateau separate heritable switching propensity from state heritability,
organ induction and founder bottlenecks?

Mechanisms simulated (all with the SAME marginal metastatic state frequency where possible):
  A  stable-state selection       eta_sd = 0, slow switching, state-dependent dissemination
  B  heritable switching          eta_sd > 0, no organ effect, no state-dependent dissemination
  C  organ induction              eta_sd = 0, organ_gain > 1 in metastatic organs
  D  bottleneck only              eta_sd = 0, no organ effect, few founders (low phi)
  E  coupled B + C

For each: the within-site, within-clone covariance curve by MRCA depth, its plateau fit, and the clone-level
variance of the state fraction. Writes results/tables/e8_*.csv.
"""
import json, os, sys, time
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.branching.simulator import Params, simulate, sample_cells, edit_paths, O
from src.lineage_inference.ctmc import mrca_depth_pairs, cov_curve, fit_plateau

QUICK = "--quick" in sys.argv
REPS = 4 if QUICK else 20
SEED = 20261124
TAB = os.path.join(ROOT, "results/tables"); MET = os.path.join(ROOT, "results/metrics")
os.makedirs(TAB, exist_ok=True); os.makedirs(MET, exist_ok=True)


def mech(name):
    """Parameter sets. Base switching a = 0.3, c = 0.7 gives a stationary high fraction of 0.3."""
    p = Params(); p.mu = 0.008          # base dissemination raised so metastases are analysable
    if name == "A_stable":
        p.a, p.b = 0.03, 0.07                       # 10x slower: state is close to heritable
        p.kappa = np.array([1.0, 8.0])              # high-EMT cells disseminate 8x more
    elif name == "B_switching":
        p.eta_sd = 0.8
    elif name == "C_induction":
        p.organ_gain = np.array([1.0, 3.0, 3.0])
    elif name == "D_bottleneck":
        p.mu = 0.001                                # 8x fewer founders (and so a smaller metastasis)
    elif name == "E_coupled":
        p.eta_sd = 0.8; p.organ_gain = np.array([1.0, 3.0, 3.0])
    else:
        raise ValueError(name)
    return p


def analyse(sim, rng, centre_within_site=True, eps=0.0):
    sm = sample_cells(sim, [4000, 4000, 4000], rng, eps=eps)
    paths = edit_paths(sm, sim["nodes"])
    rows, clone_rows = [], []
    for o in range(O):
        m = sm["organ"] == o
        if m.sum() < 120:
            continue
        st = sm["state"][m].astype(float)
        cl = sm["clone"][m]; pt = [paths[i] for i in np.where(m)[0]]
        # within-clone pairs only, centred within site (or globally)
        centre = st.mean() if centre_within_site else sm["state"].mean()
        for c in np.unique(cl):
            k = cl == c
            if k.sum() < 20:
                continue
            # centre by the site mean, not the clone mean: between-clone spread must survive in the plateau
            idx = np.where(k)[0]
            s = st[k] - centre
            dd, pp = [], []
            n = len(idx)
            pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
            if len(pairs) > 40_000:
                pairs = [pairs[t] for t in rng.choice(len(pairs), 40_000, replace=False)]
            for i, j in pairs:
                a, b = pt[idx[i]], pt[idx[j]]
                sh = 0
                for x, y in zip(a, b):
                    if x != y:
                        break
                    sh += 1
                dd.append(sh); pp.append(s[i] * s[j])
            rows += [dict(organ=o, clone=int(c), depth=h, prod=v) for h, v in zip(dd, pp)]
            clone_rows.append(dict(organ=o, clone=int(c), n=int(k.sum()), frac=float(st[k].mean())))
    return pd.DataFrame(rows), pd.DataFrame(clone_rows)


def main():
    rng = np.random.default_rng(SEED)
    summ, curves = [], []
    for name in ("A_stable", "B_switching", "C_induction", "D_bottleneck", "E_coupled"):
        for rep in range(REPS):
            p = mech(name)
            sim = simulate(p, rng)
            nmet = int((sim["organ"] > 0).sum())
            if nmet < 150:
                summ.append(dict(mech=name, rep=rep, skipped=1, n_met=nmet)); continue
            pr, cr = analyse(sim, rng)
            for site_centre in (True,):
                g = pr.groupby("depth")["prod"].agg(["mean", "size", "std"]).reset_index()
                g = g[g["size"] >= 40]
                cur = np.column_stack([g.depth, g["mean"], g["size"], g["std"] / np.sqrt(g["size"])])
                fp = fit_plateau(cur) if len(cur) >= 3 else None
                for o in sorted(pr.organ.unique()):
                    po = pr[pr.organ == o]
                    go = po.groupby("depth")["prod"].agg(["mean", "size", "std"]).reset_index()
                    go = go[go["size"] >= 40]
                    if len(go) >= 3:
                        c2 = np.column_stack([go.depth, go["mean"], go["size"], go["std"] / np.sqrt(go["size"])])
                        f2 = fit_plateau(c2)
                        curves += [dict(mech=name, rep=rep, organ=int(o), depth=int(h), cov=float(v), n=int(nn))
                                   for h, v, nn in zip(go.depth, go["mean"], go["size"])]
                    else:
                        f2 = None
                    cro = cr[cr.organ == o]
                    summ.append(dict(mech=name, rep=rep, organ=int(o), skipped=0, n_met=nmet,
                                     eta_sd=p.eta_sd, n_clones=len(cro),
                                     clone_frac_var=float(cro.frac.var(ddof=1)) if len(cro) > 1 else np.nan,
                                     mean_frac=float(cro.frac.mean()) if len(cro) else np.nan,
                                     n_founders=len(sim["founders"]),
                                     **({f"fit_{k}": v for k, v in f2.items()} if f2 else {})))
            print(name, rep, "met cells", nmet, "founders", len(sim["founders"]), flush=True)
    pd.DataFrame(summ).to_csv(os.path.join(TAB, "e8_mechanism_signatures.csv"), index=False)
    pd.DataFrame(curves).to_csv(os.path.join(TAB, "e8_cov_curves.csv"), index=False)
    print("done")


if __name__ == "__main__":
    t0 = time.time(); main(); print(f"{time.time() - t0:.0f}s")
