"""E14: interval and test calibration, with Monte Carlo uncertainty on every error rate.

Milestone 3 reported "0% false positives, 100% power" from 10 replicates and an asymptotic chi-square reference.
That reference is not trustworthy here: under H0 for a shared rate, the clone-specific rates are unidentified
when k lies outside the window, and variance components sit on a boundary. This script

  1. builds the PARAMETRIC BOOTSTRAP null distribution of each likelihood-ratio statistic by simulating under
     H0, and compares it with chi2(C-1);
  2. reports the realised type-I error of the chi-square test and of the calibrated test, each with a Wilson
     interval;
  3. measures coverage of the per-clone profile intervals for pi_c;
  4. measures coverage of a clone-bootstrap percentile interval for V_pi;
  5. repeats under a MISSPECIFIED generative process (three states collapsed to two) to test robustness.

Scenarios (C clones, balanced trees):
  H0_pi    all clones share pi; rates vary          -> the shared-pi test should reject 5% of the time
  H0_k     occupancies vary; all clones share k     -> the shared-k test should reject 5% of the time
  H0_k_frozen  as H0_k but k = 0.01 (below the identifiability window): the hard boundary case
  H1_pi / H1_k  the corresponding alternatives, for power at the calibrated threshold
  MIS      three-state generator collapsed to two states: out-of-family misspecification
"""
import os, sys, time
from multiprocessing import Pool
import numpy as np, pandas as pd
from scipy.stats import chi2

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.lineage_inference.ctmc import FastTree
from src.lineage_inference.perclone import fit_clone, variance_components, lrt_shared

TAB = os.path.join(ROOT, "results/tables"); os.makedirs(TAB, exist_ok=True)
QUICK = "--quick" in sys.argv
REPS = 40 if QUICK else 300
C = 12
DEPTH = 6
NLEAF = 64


def balanced(depth):
    par = [-1]; level = [0]
    for _ in range(depth):
        nxt = []
        for v in level:
            for _ in range(2):
                par.append(v); nxt.append(len(par) - 1)
        level = nxt
    return np.array(par), np.array(level)


def sim_tree(pi, k, rng, depth=DEPTH, three_state=False):
    """Two-state CTMC down a balanced tree. three_state collapses a 3-state chain onto {0,1} (misspecification)."""
    par, leaves = balanced(depth)
    n = len(par)
    if not three_state:
        E = np.exp(-k)
        st = np.empty(n, int); st[0] = int(rng.random() < pi)
        for i in range(1, n):
            p = st[par[i]]
            up = pi * (1 - E) if p == 0 else pi + (1 - pi) * E
            st[i] = int(rng.random() < up)
        obs = st
    else:
        # 3 states; 0 and 1 both map to observed 0, state 2 maps to observed 1. Rates chosen so the observed
        # marginal is approximately pi, but the observed process is NOT Markov.
        Pm = np.array([[1 - k * .6, k * .5, k * .1],
                       [k * .5, 1 - k * .9, k * .4],
                       [k * .15, k * .35, 1 - k * .5]])
        Pm = np.clip(Pm, 1e-6, 1); Pm /= Pm.sum(1, keepdims=True)
        st = np.empty(n, int); st[0] = rng.choice(3, p=[.4, .3, .3])
        for i in range(1, n):
            st[i] = rng.choice(3, p=Pm[st[par[i]]])
        obs = (st == 2).astype(int)
    ls = np.full(n, -1); ls[leaves] = obs[leaves]
    return FastTree(par, ls)


def make_set(scen, rng):
    trees, truth = [], []
    for c in range(C):
        if scen == "H0_pi":
            pi, k = 0.5, float(np.exp(rng.normal(np.log(0.5), 0.9)))
        elif scen == "H1_pi":
            pi, k = float(np.clip(rng.normal(0.5, 0.18), .05, .95)), 0.5
        elif scen == "H0_k":
            pi, k = float(np.clip(rng.normal(0.5, 0.18), .05, .95)), 0.5
        elif scen == "H0_k_frozen":
            pi, k = float(np.clip(rng.normal(0.5, 0.18), .05, .95)), 0.01
        elif scen == "H1_k":
            pi, k = 0.5, float(np.exp(rng.normal(np.log(0.5), 0.9)))
        elif scen == "MIS":
            pi, k = 0.5, 0.5
        else:
            raise ValueError(scen)
        trees.append(sim_tree(pi, k, rng, three_state=(scen == "MIS")))
        truth.append((pi, k))
    return trees, np.array(truth)


def one(args):
    scen, rep = args
    rng = np.random.default_rng(770000 + 7919 * rep + abs(hash(scen)) % 997)
    trees, truth = make_set(scen, rng)
    fits = [fit_clone(t) for t in trees]
    vc = variance_components(fits, min_n=10)
    lp = lrt_shared(trees, fits, "pi"); lk = lrt_shared(trees, fits, "k")
    pis = np.array([f["pi"] for f in fits])
    cov_pi = float(np.mean([(f["pi_lo"] - 1e-9) <= t <= (f["pi_hi"] + 1e-9)
                            for f, t in zip(fits, truth[:, 0])]))
    # clone-bootstrap percentile interval for V_pi, resampling the per-clone estimates
    se = np.array([f["pi_se"] if np.isfinite(f["pi_se"]) else 0.0 for f in fits])
    bs = []
    for _ in range(400):
        idx = rng.choice(len(pis), len(pis), replace=True)
        bs.append(max(np.var(pis[idx], ddof=1) - np.mean(se[idx] ** 2), 0.0))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    tv_pi = float(np.var(truth[:, 0], ddof=1)); tv_k = float(np.var(np.log(truth[:, 1]), ddof=1))
    return dict(scen=scen, rep=rep, lrt_pi=lp["stat"], lrt_k=lk["stat"],
                p_chi2_pi=lp["p_asym"], p_chi2_k=lk["p_asym"],
                true_V_pi=tv_pi, true_V_logk=tv_k, V_pi=vc.get("V_pi", np.nan),
                V_pi_lo=float(lo), V_pi_hi=float(hi),
                cover_V_pi=float(lo <= tv_pi <= hi), cover_pi_profile=cov_pi,
                n_k_ident=vc.get("n_clones_k_identified", 0))


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def main():
    scens = ["H0_pi", "H1_pi", "H0_k", "H0_k_frozen", "H1_k", "MIS"]
    jobs = [(s, r) for s in scens for r in range(REPS)]
    with Pool(12) as p:
        rows = p.map(one, jobs)
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(TAB, "e14_calibration_raw.csv"), index=False)

    out = []
    # calibrated thresholds from the parametric-bootstrap nulls
    thr_pi = np.percentile(d[d.scen == "H0_pi"].lrt_pi, 95)
    thr_k = np.percentile(d[d.scen == "H0_k"].lrt_k, 95)
    thr_k_frozen = np.percentile(d[d.scen == "H0_k_frozen"].lrt_k, 95)
    for s, g in d.groupby("scen"):
        n = len(g)
        for name, rej in [("chi2_pi", (g.p_chi2_pi < .05)), ("chi2_k", (g.p_chi2_k < .05)),
                          ("calib_pi", (g.lrt_pi > thr_pi)), ("calib_k", (g.lrt_k > thr_k))]:
            kk = int(rej.sum()); lo, hi = wilson(kk, n)
            out.append(dict(scen=s, test=name, n=n, rate=kk / n, lo=lo, hi=hi))
        out.append(dict(scen=s, test="cover_V_pi", n=n, rate=g.cover_V_pi.mean(),
                        **dict(zip(["lo", "hi"], wilson(int(g.cover_V_pi.sum()), n)))))
        out.append(dict(scen=s, test="cover_pi_profile", n=n, rate=g.cover_pi_profile.mean(),
                        lo=np.nan, hi=np.nan))
    o = pd.DataFrame(out)
    o.to_csv(os.path.join(TAB, "e14_calibration_summary.csv"), index=False)
    print(o.round(4).to_string())
    print(f"\ncalibrated 95% thresholds: pi {thr_pi:.2f}, k {thr_k:.2f}, k(frozen null) {thr_k_frozen:.2f}; "
          f"chi2(df={C - 1}) 95% = {chi2.ppf(.95, C - 1):.2f}")
    q = pd.DataFrame(dict(scen=["H0_pi", "H0_k", "H0_k_frozen"],
                          thr_boot=[thr_pi, thr_k, thr_k_frozen],
                          thr_chi2=[chi2.ppf(.95, C - 1)] * 3))
    q.to_csv(os.path.join(TAB, "e14_thresholds.csv"), index=False)


if __name__ == "__main__":
    t0 = time.time(); main(); print(f"{time.time() - t0:.0f}s")
