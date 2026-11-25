"""E15b: recorder-capacity design law, and the occupancy-versus-kinetics budget trade-off.

The first formulation (run_e15.py) used the Fisher information about the grand occupancy, which is fixed by the
cell budget and so could not discriminate designs, and modelled staggering in a way that could only lose
resolution. Both are corrected here.

Design model
------------
A clone is grown for D cell divisions and the recorder has S usable recording sites, so the reconstructed tree
resolves D_eff = min(S, D) coalescent levels, each spanning b = D / D_eff divisions. Lineage distances available
to the analyst are therefore L in {b, 2b, ..., D}.

From the information calculus (E12), a pair at distance L is informative about k when k*L lies in roughly
[0.10, 1.0]. The set of rates some pair can speak to is therefore

        k in [ 0.10 / D ,  1.0 / b ] = [ 0.10 / D , D_eff / D ],

so the number of decades of switching rate a recorder covers is

        W = log10( 10 * D_eff ) = 1 + log10(D_eff).

The upper limit is set by the FINEST resolvable branch (one edit), not by tree depth; the lower limit by the
total recording duration. Objectives:
  * I_Vpi    information about the occupancy variance, which scales with the number of clones;
  * I_logk   worst-case and geometric-mean information about the rate over an uncertainty range;
  * W        decades of rate covered.
"""
import os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from experiments.sampling_design.fisher_audit import fisher2

TAB = os.path.join(ROOT, "results/tables"); os.makedirs(TAB, exist_ok=True)
BUDGET = 2048
K_RANGE = np.exp(np.linspace(np.log(0.02), np.log(5.0), 31))


def pair_info(D_eff, b, k, pi):
    """Total pair Fisher information about (log k, pi) for one balanced resolved tree."""
    I = np.zeros((2, 2))
    for h in range(D_eff):
        npair = 2 ** h * (2 ** (D_eff - h - 1)) ** 2
        I += npair * fisher2(pi, k, (D_eff - h) * b)
    return I


def info_Vpi(C, n, pi, k, b, D_eff):
    """Information about Var_c(pi_c): a variance over C clones, each measured with finite precision.

    Var(sample variance) ~ 2 (V + s^2)^2 / (C - 1) where s^2 is the per-clone sampling variance of pi_hat.
    Within-clone correlation inflates s^2 by the design effect 1 + (m - 1) rho_bar with m the mean cluster size.
    """
    if C < 3:
        return 0.0
    rho_bar = np.mean([np.exp(-2 * k * (D_eff - h) * b) for h in range(D_eff)])
    deff = 1 + (n - 1) * max(rho_bar, 0.0)
    s2 = pi * (1 - pi) / max(n / deff, 1.0)
    V = 0.03                                     # a representative occupancy variance
    return (C - 1) / (2 * (V + s2) ** 2)


def main():
    rows = []
    for D in (6, 8, 10, 12):
        for S in (2, 4, 8, 16, 32):
            D_eff = min(S, D); b = D / D_eff
            for n in (32, 64, 128, 256, 512):
                C = BUDGET // n
                if C < 3:
                    continue
                for pi in (0.5,):
                    ik, ivp = [], []
                    for k in K_RANGE:
                        I = pair_info(D_eff, b, k, pi) * C
                        ik.append(I[0, 0] - (I[0, 1] ** 2 / I[1, 1] if I[1, 1] > 1e-12 else 0.0))
                        ivp.append(info_Vpi(C, n, pi, k, b, D_eff))
                    ik = np.array(ik)
                    rows.append(dict(D=D, S=S, D_eff=D_eff, branch=b, cells_per_clone=n, C=C, pi=pi,
                                     k_lo_cov=0.10 / D, k_hi_cov=D_eff / D,
                                     decades_covered=np.log10(max(D_eff, 1e-9) / 0.10),
                                     I_logk_worst=ik.min(),
                                     I_logk_geomean=float(np.exp(np.mean(np.log(np.clip(ik, 1e-12, None))))),
                                     I_Vpi=float(np.median(ivp))))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(TAB, "e15b_capacity_designs.csv"), index=False)

    print("decades of switching rate covered, by recorder capacity (independent of cell budget):")
    print(d.groupby(["S", "D"]).decades_covered.first().unstack().round(2).to_string())

    print("\nbest design per objective (budget %d cells):" % BUDGET)
    for obj in ("I_logk_worst", "I_logk_geomean", "I_Vpi"):
        b_ = d.loc[d[obj].idxmax()]
        print(f"  {obj:16s} -> S={int(b_.S)} sites, D={int(b_.D)} divisions, {int(b_.C)} clones x "
              f"{int(b_.cells_per_clone)} cells  (value {b_[obj]:.3g})")

    # Pareto frontier between occupancy-variance and worst-case rate information
    par = []
    for _, r in d.iterrows():
        dom = ((d.I_Vpi >= r.I_Vpi) & (d.I_logk_worst >= r.I_logk_worst) &
               ((d.I_Vpi > r.I_Vpi) | (d.I_logk_worst > r.I_logk_worst))).any()
        par.append(not dom)
    d["pareto"] = par
    d.to_csv(os.path.join(TAB, "e15b_capacity_designs.csv"), index=False)
    pf = d[d.pareto].sort_values("I_Vpi")
    print(f"\nPareto-optimal designs: {len(pf)} of {len(d)}")
    print(pf[["S", "D", "C", "cells_per_clone", "decades_covered", "I_logk_worst", "I_Vpi"]]
          .round(3).to_string(index=False))

    # where do the published recorders sit?
    real = pd.DataFrame([dict(system="macsGESTALT (Simeonov)", D_eff=4, D=10),
                         dict(system="Cassiopeia (Quinn)", D_eff=12, D=14)])
    real["k_hi_cov"] = real.D_eff / real.D
    real["k_lo_cov"] = 0.10 / real.D
    real["decades"] = np.log10(real.D_eff / 0.10)
    real.to_csv(os.path.join(TAB, "e15b_real_recorders.csv"), index=False)
    print("\npublished recorders, on the same scale:")
    print(real.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
