"""Table B6: what each method establishes, and what it leaves unresolved.

Qualitative, but every cell is traceable to an executed result; the `evidence` column
names the table it comes from.

Usage: python -I make_tableB6.py <benchdir>
"""
import os, sys, csv
B = sys.argv[1]
TAB = os.path.join(B, "results", "tables")
os.makedirs(TAB, exist_ok=True)

rows = [
 dict(method="Study D", establishes="per-clone occupancy and its across-clone variance, with "
      "calibrated tests and profile intervals; an explicit identifiability window for the rate",
      unresolved="over-disperses per-clone rates: reports the largest spurious V_log k of the "
      "three binary methods when the true value is zero",
      unique_contribution="states which rates the tree can resolve, and flags its own boundary cases",
      evidence="TableB3, TableB3b, TableB4b"),
 dict(method="PATH", establishes="phylogenetic autocorrelation with an analytical z-test; a "
      "transition matrix Q whose stationary distribution recovers occupancy",
      unresolved="rate estimate is nearly flat in the true rate and biased low by 1.6-2.3x; "
      "no interval on Q; fails on clones fixed in one state",
      unique_contribution="the only method here with an analytical significance test for "
      "phylogenetic signal itself",
      evidence="TableB3, TableB4b, Figure B3C"),
 dict(method="PASTRI", establishes="depth-stratified transition matrices, so rates are allowed "
      "to change with developmental stage; occupancy recovered as well as the others",
      unresolved="rate biased high by 11-16x and exactly invariant to branch times; 3-18% of "
      "clones return no finite rate; no per-rate interval",
      unique_contribution="the only method that does not assume a single rate across the tree",
      evidence="TableB3, TableB4, Figure B3A"),
 dict(method="SCOUT", establishes="OU parameters for a continuous trait with AICc model "
      "selection; theta tracks clone occupancy in most mechanisms",
      unresolved="theta carries no signal when states are near-frozen; alpha is in units of "
      "whatever clock the tree carries, which is unit edges when the tree has none",
      unique_contribution="the only method that models a continuous phenotype natively, and "
      "the only one that selects between evolutionary models",
      evidence="TableB3, method_reproduction.csv"),
 dict(method="All four, jointly", establishes="clone-associated occupancy is a real and "
      "recoverable property of lineage-recorded data, and the finding does not depend on which "
      "method is used",
      unresolved="switching kinetics: three implementations disagree by more than an order of "
      "magnitude on identical trees, bracketing the truth from both sides",
      unique_contribution="the disagreement itself is the evidence that the rate is not "
      "identified by these data",
      evidence="TableB3, Figure B2B"),
]
f = os.path.join(TAB, "TableB6_scientific_contribution.csv")
with open(f, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print("wrote", f)
for r in rows:
    print(f"\n### {r['method']}\n  establishes : {r['establishes']}\n  unresolved  : {r['unresolved']}")
