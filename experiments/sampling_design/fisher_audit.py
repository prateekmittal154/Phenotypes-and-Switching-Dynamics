"""E12: audit of the Milestone 3 Fisher-information claim kL* ~ 0.405.

The Milestone 3 statement was imprecise on four points the Milestone 4 brief raises. This script fixes the
definitions and recomputes each case.

Definitions used here
---------------------
* Parameter of interest: log k, where k = a + b is the total switching rate of the two-state CTMC.
* Observation: the joint state (X_i, X_j) of one within-clone leaf pair, a single draw from a 4-cell multinomial.
* L: the lineage distance from the MRCA to EACH leaf, so the total path between the pair is 2L and the state
  correlation is rho = exp(-2kL). Milestone 3 used this convention but did not state it; in terms of the total
  pair distance d = 2L the optimum is k*d ~ 0.81.
* pi known vs estimated: Milestone 3 computed the information about log k with pi held FIXED. When pi is
  estimated jointly, the relevant quantity is the Schur complement I_kk - I_k,pi^2 / I_pi,pi.
* Independence: Milestone 3 treated pairs as independent, which is false inside a tree. Quantified in E13.

Outputs results/tables/e12_*.csv
"""
import os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
TAB = os.path.join(ROOT, "results/tables"); os.makedirs(TAB, exist_ok=True)


def probs(pi, k, L):
    q = 1 - pi
    rho = np.exp(-2 * k * L)
    return np.array([q * q + pi * q * rho, pi * q * (1 - rho), pi * q * (1 - rho), pi * pi + pi * q * rho])


def fisher2(pi, k, L, h=1e-6):
    """Full 2x2 Fisher information for (log k, pi) from one pair, by numerical differentiation."""
    p = np.clip(probs(pi, k, L), 1e-14, None)
    dk = (probs(pi, k * np.exp(h), L) - probs(pi, k * np.exp(-h), L)) / (2 * h)
    dp = (probs(pi + h, k, L) - probs(pi - h, k, L)) / (2 * h)
    I = np.zeros((2, 2))
    I[0, 0] = (dk * dk / p).sum(); I[1, 1] = (dp * dp / p).sum(); I[0, 1] = I[1, 0] = (dk * dp / p).sum()
    return I


def main():
    rows = []
    for k in (0.05, 0.1, 0.5, 2.0, 5.0):
        for pi in (0.1, 0.3, 0.5, 0.7):
            for L in np.exp(np.linspace(np.log(1e-3), np.log(200), 400)):
                I = fisher2(pi, k, L)
                known = I[0, 0]
                joint = I[0, 0] - (I[0, 1] ** 2 / I[1, 1] if I[1, 1] > 1e-12 else 0.0)
                rows.append(dict(k=k, pi=pi, L=L, kL=k * L, I_logk_pi_known=known, I_logk_pi_joint=joint))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(TAB, "e12_fisher_full.csv"), index=False)
    opt = []
    for (k, pi), g in d.groupby(["k", "pi"]):
        a = g.loc[g.I_logk_pi_known.idxmax()]
        b = g.loc[g.I_logk_pi_joint.idxmax()]
        opt.append(dict(k=k, pi=pi, kL_opt_pi_known=a.kL, kL_opt_pi_joint=b.kL,
                        I_max_known=a.I_logk_pi_known, I_max_joint=b.I_logk_pi_joint,
                        ratio_joint_over_known=b.I_logk_pi_joint / a.I_logk_pi_known))
    o = pd.DataFrame(opt)
    o.to_csv(os.path.join(TAB, "e12_optimum_audit.csv"), index=False)
    print(o.round(4).to_string(), flush=True)
    print("\nkL* (pi known):  mean %.4f  sd %.4f  range %.4f-%.4f"
          % (o.kL_opt_pi_known.mean(), o.kL_opt_pi_known.std(), o.kL_opt_pi_known.min(), o.kL_opt_pi_known.max()))
    print("kL* (pi joint):  mean %.4f  sd %.4f  range %.4f-%.4f"
          % (o.kL_opt_pi_joint.mean(), o.kL_opt_pi_joint.std(), o.kL_opt_pi_joint.min(), o.kL_opt_pi_joint.max()))

    # how flat is the optimum? the practical question for a designer
    flat = []
    for (k, pi), g in d.groupby(["k", "pi"]):
        m = g.I_logk_pi_joint.max()
        ok = g[g.I_logk_pi_joint >= 0.5 * m]
        flat.append(dict(k=k, pi=pi, kL_lo_50pct=ok.kL.min(), kL_hi_50pct=ok.kL.max(),
                         fold_width=ok.kL.max() / max(ok.kL.min(), 1e-9)))
    f = pd.DataFrame(flat); f.to_csv(os.path.join(TAB, "e12_optimum_width.csv"), index=False)
    print("\nrange of k*L retaining >=50% of peak information:")
    print(f.round(3).to_string())


if __name__ == "__main__":
    main()
