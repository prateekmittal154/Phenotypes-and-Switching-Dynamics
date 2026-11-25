"""E15: budgeted design optimisation, including robust (minimax) designs when the switching rate is unknown.

Design model
------------
A clone is grown for D cell divisions, giving n = 2^D sampled cells. The recorder marks a division with
probability p, so the reconstructed tree resolves D_eff = pD coalescent levels and one resolved branch spans
b = 1/p divisions. Lineage distance between a pair whose resolved MRCA sits at level h is L = (D_eff - h) * b,
measured in divisions. Lowering p therefore buys longer resolvable distances at coarser granularity, which is
exactly the knob a recorder designer controls (edit rate, induction window, barcode capacity).

Total information for a design is the sum over within-clone pair classes of (number of pairs) x (Fisher
information of a pair at that distance), using the per-pair information of Proposition U5. Pairs in a tree are
not independent, so this overstates the true information; E13 measures the inflation factor and it is applied
here as a constant deflation, which is conservative for comparing designs at equal budget.

Objectives compared: information about log k; about pi; the worst case of the former over an uncertainty range
for k; and the Pareto frontier between the two.
"""
import itertools, os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from experiments.sampling_design.fisher_audit import fisher2

TAB = os.path.join(ROOT, "results/tables"); os.makedirs(TAB, exist_ok=True)
BUDGET = 2048                      # total sampled cells
K_RANGE = np.exp(np.linspace(np.log(0.02), np.log(3.0), 25))   # plausible rates, unknown to the designer


def pair_classes(D_eff, b):
    """(lineage distance, number of pairs) for a balanced resolved tree of depth D_eff with branch length b."""
    D_eff = int(round(D_eff))
    if D_eff < 1:
        return []
    out = []
    n = 2 ** D_eff
    for h in range(0, D_eff):
        # pairs whose MRCA is at level h: each of 2^h nodes contributes (2^(D_eff-h-1))^2 cross pairs
        npair = 2 ** h * (2 ** (D_eff - h - 1)) ** 2
        L = (D_eff - h) * b
        out.append((L, npair))
    return out


def design_info(C, D, p, k, pi, deflate=1.0):
    """Total pair information about (log k, pi) for C clones, each 2^D cells, recorder edit probability p."""
    D_eff = max(1, int(round(p * D)))
    b = 1.0 / p
    I = np.zeros((2, 2)); npairs = 0
    for L, m in pair_classes(D_eff, b):
        I += m * fisher2(pi, k, L); npairs += m
    I *= C / max(deflate, 1e-9)
    joint_k = I[0, 0] - (I[0, 1] ** 2 / I[1, 1] if I[1, 1] > 1e-12 else 0.0)
    # information about pi also accrues from the leaf marginal across clones (C * n Bernoulli draws)
    n = 2 ** D
    I_pi_marginal = C * n / max(pi * (1 - pi), 1e-9)
    return joint_k, I_pi_marginal, npairs * C, D_eff, b


def main():
    deflate = 1.0
    f13 = os.path.join(TAB, "e13_tree_vs_pair_mc.csv")
    if os.path.exists(f13):
        deflate = float(pd.read_csv(f13).inflation.median())
        print(f"applying naive-pair deflation factor {deflate:.1f}x from E13", flush=True)

    rows = []
    for D in (5, 6, 7, 8, 9):
        n = 2 ** D
        C = BUDGET // n
        if C < 3:
            continue
        for p in (0.1, 0.2, 0.35, 0.5, 0.75, 1.0):
            for pi in (0.3, 0.5):
                vals = []
                for k in K_RANGE:
                    jk, ip, npair, De, b = design_info(C, D, p, k, pi, deflate)
                    vals.append((k, jk, ip))
                v = np.array(vals)
                rows.append(dict(D=D, cells_per_clone=n, C=C, p_edit=p, pi=pi, D_eff=De, branch=b,
                                 I_logk_worst=v[:, 1].min(), I_logk_geomean=float(np.exp(np.mean(np.log(
                                     np.clip(v[:, 1], 1e-12, None))))),
                                 I_logk_at_k0p5=float(np.interp(0.5, v[:, 0], v[:, 1])),
                                 I_pi=v[0, 2], n_pairs=npair))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(TAB, "e15_designs.csv"), index=False)

    print("\nbest designs by objective (pi = 0.5):")
    g = d[d.pi == 0.5]
    for obj in ("I_logk_worst", "I_logk_geomean", "I_logk_at_k0p5", "I_pi"):
        b = g.loc[g[obj].idxmax()]
        print(f"  {obj:18s} -> {int(b.C)} clones x {int(b.cells_per_clone)} cells, p_edit={b.p_edit}, "
              f"resolved depth {int(b.D_eff)}, branch {b.branch:.1f}  value {b[obj]:.3g}")

    # staggered recording: split the budget across several edit rates
    stag = []
    for combo in itertools.combinations((0.1, 0.2, 0.35, 0.5, 1.0), 2):
        for D in (6, 7, 8):
            n = 2 ** D; C = BUDGET // n
            if C < 4:
                continue
            half = C // 2
            worst = np.inf; vals = []
            for k in K_RANGE:
                tot = sum(design_info(half, D, p, k, 0.5, deflate)[0] for p in combo)
                vals.append(tot); worst = min(worst, tot)
            stag.append(dict(mix=f"{combo[0]}+{combo[1]}", D=D, C=C,
                             I_logk_worst=worst, I_logk_geomean=float(np.exp(np.mean(np.log(
                                 np.clip(vals, 1e-12, None)))))))
    s = pd.DataFrame(stag)
    s.to_csv(os.path.join(TAB, "e15_staggered.csv"), index=False)
    best_single = g.I_logk_worst.max()
    print(f"\nbest single-rate worst-case information:  {best_single:.3g}")
    if len(s):
        bs = s.loc[s.I_logk_worst.idxmax()]
        print(f"best staggered (two edit rates):          {bs.I_logk_worst:.3g}  "
              f"({bs['mix']}, depth {int(bs.D)}, {int(bs.C)} clones)")
        print(f"gain from staggering: {bs.I_logk_worst / max(best_single, 1e-12):.2f}x worst-case")

    # Pareto frontier between occupancy and switching information
    par = []
    for _, r in g.iterrows():
        dominated = ((g.I_logk_worst >= r.I_logk_worst) & (g.I_pi >= r.I_pi) &
                     ((g.I_logk_worst > r.I_logk_worst) | (g.I_pi > r.I_pi))).any()
        par.append(not dominated)
    g2 = g.copy(); g2["pareto"] = par
    g2.to_csv(os.path.join(TAB, "e15_pareto.csv"), index=False)
    print(f"\nPareto-optimal designs: {int(g2.pareto.sum())} of {len(g2)}")
    print(g2[g2.pareto][["C", "cells_per_clone", "p_edit", "D_eff", "I_logk_worst", "I_pi"]]
          .sort_values("I_pi").round(3).to_string())


if __name__ == "__main__":
    main()
