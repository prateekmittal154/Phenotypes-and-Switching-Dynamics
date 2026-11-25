"""E9: mechanism recovery with the two-step decomposition, plus misspecification and nuisance stress tests.

For each mechanism (A stable state, B heritable propensity, C organ induction, D founder bottleneck, E coupled)
the estimator is applied to the primary tumour and to each metastatic site. Ground truth for B is
Var_c(pi_c), computed exactly from the simulated clone propensities.

Stress tests, all on mechanism B: state misclassification eps, lost clone labels, and fewer sampled cells.
"""
import os, sys, time
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.branching.simulator import Params, simulate, sample_cells, edit_paths, O
from src.lineage_inference.propensity import decompose, curve_frame

QUICK = "--quick" in sys.argv
REPS = 3 if QUICK else 15
NBOOT = 40 if QUICK else 150
SEED = 20261125
TAB = os.path.join(ROOT, "results/tables"); os.makedirs(TAB, exist_ok=True)


def mech(name):
    p = Params(); p.mu = 0.008
    if name == "A_stable":
        p.a, p.b = 0.03, 0.07; p.kappa = np.array([1.0, 8.0])
    elif name == "B_switching":
        p.eta_sd = 0.8
    elif name == "C_induction":
        p.organ_gain = np.array([1.0, 3.0, 3.0])
    elif name == "D_bottleneck":
        p.mu = 0.001
    elif name == "E_coupled":
        p.eta_sd = 0.8; p.organ_gain = np.array([1.0, 3.0, 3.0])
    return p


def truth_B(p):
    """Var over clones of the stationary high-state fraction pi_c = a_c / (a_c + b_c) in the primary."""
    if p.eta is None or p.eta_sd == 0:
        return 0.0
    a = p.a * np.exp(p.eta); b = p.b * np.exp(-p.eta)
    return float(np.var(a / (a + b), ddof=1))


def run(name, rep, rng, eps=0.0, clone_loss=0.0, n_cells=4000, nboot=NBOOT):
    p = mech(name)
    sim = simulate(p, rng)
    sm = sample_cells(sim, [n_cells] * O, rng, eps=eps, clone_loss=clone_loss)
    paths = edit_paths(sm, sim["nodes"])
    rows = []
    for o in range(O):
        m = sm["organ"] == o
        if m.sum() < 120:
            continue
        sub = [paths[i] for i in np.where(m)[0]]
        r = decompose(sub, sm["state"][m], sm["clone"][m], n_boot=nboot, rng=rng)
        if r:
            rows.append(dict(mech=name, rep=rep, organ=o, eps=eps, clone_loss=clone_loss, n_sample=n_cells,
                             truth_B=truth_B(p), n_founders=len(sim["founders"]), **r))
    return rows, sim, sm, paths


def main():
    rng = np.random.default_rng(SEED)
    rows, cur = [], []
    for name in ("A_stable", "B_switching", "C_induction", "D_bottleneck", "E_coupled"):
        for rep in range(REPS):
            r, sim, sm, paths = run(name, rep, rng)
            rows += r
            if rep == 0:
                for o in range(O):
                    m = sm["organ"] == o
                    if m.sum() < 120:
                        continue
                    cf = curve_frame([paths[i] for i in np.where(m)[0]], sm["state"][m], sm["clone"][m], rng=rng)
                    if cf is not None:
                        cur += [dict(mech=name, organ=o, depth=int(a), clone_centred=float(b),
                                     site_centred=float(c), n_pairs=int(d)) for a, b, c, d in cf]
            print(name, rep, flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(TAB, "e9_mechanism_recovery.csv"), index=False)
    pd.DataFrame(cur).to_csv(os.path.join(TAB, "e9_curves.csv"), index=False)

    stress = []
    for eps in (0.0, 0.05, 0.15):
        for cl in (0.0, 0.3):
            for nc in (4000, 800):
                for rep in range(max(3, REPS // 3)):
                    r, *_ = run("B_switching", rep, rng, eps=eps, clone_loss=cl, n_cells=nc, nboot=40)
                    stress += r
                print("stress", eps, cl, nc, flush=True)
    pd.DataFrame(stress).to_csv(os.path.join(TAB, "e9_stress.csv"), index=False)
    print("done")


if __name__ == "__main__":
    t0 = time.time(); main(); print(f"{time.time() - t0:.0f}s")
