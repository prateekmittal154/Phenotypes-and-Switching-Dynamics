"""E10-focus: the single decisive comparison, at one clone count, run in parallel over mechanisms.

OCC and FROZE share an occupancy distribution and differ only in switching timescale. The Milestone 2 plateau
cannot separate them by construction (U1). This asks whether the per-clone estimator can.
"""
import os, sys, time
from multiprocessing import Pool
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.lineage_inference.perclone import fit_clone, variance_components, lrt_shared
from experiments.lineage_identifiability.run_e10 import make

REPS = 10
C = 20


def one(args):
    mech, rep = args
    rng = np.random.default_rng(900000 + 7919 * rep + hash(mech) % 1000)
    trees, truth = make(mech, C, rng)
    fits = [fit_clone(t) for t in trees]
    vc = variance_components(fits)
    lp = lrt_shared(trees, fits, "pi"); lk = lrt_shared(trees, fits, "k")
    return dict(mech=mech, rep=rep, C=C,
                true_V_pi=float(np.var(truth[:, 0], ddof=1)),
                true_V_logk=float(np.var(np.log(truth[:, 1]), ddof=1)),
                lrt_pi_p=lp["p_asym"], lrt_k_p=lk["p_asym"], **vc)


if __name__ == "__main__":
    t0 = time.time()
    jobs = [(m, r) for m in ("OCC", "RATE", "FROZE", "BOTH") for r in range(REPS)]
    with Pool(12) as p:
        rows = p.map(one, jobs)
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(ROOT, "results/tables/e10_focus.csv"), index=False)
    cols = [c for c in ["true_V_pi", "V_pi", "true_V_logk", "V_logk", "n_clones_k_identified",
                        "lrt_pi_p", "lrt_k_p"] if c in d.columns]
    print(d.groupby("mech")[cols].median().round(4).to_string())
    print(d.assign(sig_pi=d.lrt_pi_p < .05, sig_k=d.lrt_k_p < .05).groupby("mech")[
        ["sig_pi", "sig_k"]].mean().round(3).to_string())
    print(f"{time.time() - t0:.0f}s")
