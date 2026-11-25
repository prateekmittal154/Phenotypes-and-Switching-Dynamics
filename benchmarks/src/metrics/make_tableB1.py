"""Table B1: the estimand and assumptions matrix, machine readable.

Mirrors protocols/estimand_matrix.md, which was fixed before any run.
Usage: python -I make_tableB1.py <benchdir>
"""
import os, sys, csv
B = sys.argv[1]
TAB = os.path.join(B, "results", "tables"); os.makedirs(TAB, exist_ok=True)

rows = [
 dict(attribute="Native estimand",
      StudyD="per-clone (pi_c, k_c) by ML; across-clone variance components V_pi, V_log k",
      PATH="phylogenetic auto/cross-correlation (xcor); transition matrix P and rate matrix Q=logm(P)",
      PASTRI="per-cell-cycle transition probability matrix T, per depth window d<=fi_depth",
      SCOUT="Ornstein-Uhlenbeck alpha (reversion rate), theta (optimum), sigma^2 (drift)"),
 dict(attribute="Input phenotype", StudyD="binary", PATH="categorical (one-hot) or continuous",
      PASTRI="categorical cell-type labels", SCOUT="continuous"),
 dict(attribute="Tree requirements",
      StudyD="topology and per-edge lengths (VarTree), or unit edges (FastTree)",
      PATH="topology; branch lengths enter only via the mean patristic distance",
      PASTRI="topology only; depth is node-count based (lca_normalized_height)",
      SCOUT="topology; branch lengths if present, unit edges substituted if absent"),
 dict(attribute="Dynamics assumed", StudyD="two-state time-homogeneous CTMC per clone",
      PATH="time-homogeneous Markov transitions; stationarity for the correlation statistic",
      PASTRI="rates may vary with depth; optional hard irreversibility constraints",
      SCOUT="OU process; regimes may be assigned per clade"),
 dict(attribute="Output type", StudyD="point estimates, profile intervals, calibrated LR tests",
      PATH="point estimates; analytical z-test for xcor only; no interval on Q",
      PASTRI="point estimate per depth window; permutation test in the paper, no per-rate CI",
      SCOUT="point estimates plus AICc model selection"),
 dict(attribute="Compatible simulations", StudyD="S1-S8 (binary)", PATH="S1-S8 (binary)",
      PASTRI="S1-S8 (binary); S5 targets its stage-specific design",
      SCOUT="the OU trait generated alongside the binary state in every scenario"),
 dict(attribute="Valid Study D comparator", StudyD="n/a",
      PATH="Q off-diagonals <-> k; stationary distribution of Q <-> pi",
      PASTRI="off-diagonals of T, converted by k=-log(1-p) <-> k; stationary distribution <-> pi",
      SCOUT="alpha <-> k ONLY as an autocorrelation decay rate; Var_c(theta) <-> V_pi"),
 dict(attribute="Comparisons that must NOT be made", StudyD="n/a",
      PATH="no interval-coverage comparison for PATH_inf; absolute rate only on a matched branch-length convention",
      PASTRI="no interval coverage; invariance to branch rescaling is a design consequence, not an error",
      SCOUT="alpha is not a binary switching rate and is never plotted on the same axis as k"),
 dict(attribute="Reports uncertainty usable for coverage", StudyD="yes (profile likelihood)",
      PATH="no (for PATH_inf)", PASTRI="no", SCOUT="yes (likelihood / AICc)"),
]
f = os.path.join(TAB, "TableB1_estimand_matrix.csv")
with open(f, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["attribute", "StudyD", "PATH", "PASTRI", "SCOUT"])
    w.writeheader(); w.writerows(rows)
print("wrote", f, f"({len(rows)} attributes)")
