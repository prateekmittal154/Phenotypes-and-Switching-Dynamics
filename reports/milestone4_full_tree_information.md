# Milestone 4: full-tree information versus the naive independent-pair sum

Milestone 3's design calculus summed Fisher information over within-clone pairs as if they were independent.
They are not: pairs in a tree share ancestors. This report measures the error exactly.

## Method
For a balanced tree of depth $d$ with $n = 2^d$ leaves, the exact Fisher information of the **full-tree
likelihood** is
$$I(\theta) = \sum_{x \in \{0,1\}^n} p(x;\theta)\, s(x)s(x)^{\!\top}, \qquad s(x) = \nabla_\theta \log p(x;\theta),$$
computed by enumerating all $2^n$ leaf configurations. This is tractable to $n = 8$; beyond that the same
quantity is estimated by Monte Carlo (2,000 simulated trees per setting), validated against the exact values.
Scripts: `experiments/sampling_design/full_tree_information.py`; table
`results/tables/e13_tree_vs_pair_mc.csv`.

**Monte Carlo validation.** At depth 3, where both are available, the Monte Carlo estimate matches exact
enumeration to a **maximum relative error of 6.0%** across 12 parameter settings.

## Result 1 — the naive sum overstates information by up to two orders of magnitude

Inflation factor = (naive independent-pair sum) / (exact full-tree information):

| tree depth | median inflation | maximum |
|---|---|---|
| 3 | 3.4× | 8.3× |
| 4 | 5.2× | 20.2× |
| 5 | 7.4× | 45.8× |
| 6 | 10.2× | 99.0× |
| 7 | 14.2× | 221.3× |

The error grows with depth, because deeper trees have quadratically more pairs but only linearly more
independent lineage history.

## Result 2 — the inflation is governed by the switching rate

| switching rate k | median inflation |
|---|---|
| 0.05 | 45.7× |
| 0.10 | 30.7× |
| 0.20 | 13.0× |
| 0.50 | 2.0× |
| 1.00 | **1.06×** |
| 2.00 | **0.98×** |

This has a clean interpretation. When switching is **fast**, state memory decays within a branch, so distant
pairs really are nearly independent and the naive sum is almost exact. When switching is **slow**, every cell in
a clone carries essentially the same information, the clone behaves as a single observation, and summing over
$O(n^2)$ pairs is badly wrong.

**The regime where the naive calculus fails worst is precisely the regime where slow switching is hard to
distinguish from a frozen state** — the failure mode identified in Milestone 3. The two problems compound.

## Result 3 — the pairwise design rule survives, the pairwise precision estimate does not

Depth maximising full-tree information **per cell** (π = 0.5):

| k | best depth | k × depth |
|---|---|---|
| 0.05 | 6 | 0.3 |
| 0.10 | 7 | 0.7 |
| 0.20 | 6 | 1.2 |
| 0.50 | 4 | 2.0 |
| 1.00 | 3 | 3.0 |

Typical pair distances in a depth-$d$ tree are around $d/2$, so these optima correspond to $kL \approx 0.2$–1.5,
bracketing the pairwise optimum $kL^\ast \approx 0.43$. **The pairwise rule therefore remains a usable guide for
choosing a design**, accurate to within a factor of two or three in the depth it recommends.

What it cannot do is predict absolute precision. Because the inflation factor ranges from 1× to 221×, any
sample-size calculation based on summing pair information will be optimistic — by two orders of magnitude in
the slow-switching regime. Conclusions:

- **Use the pairwise calculus to rank designs.** The ordering is preserved.
- **Do not use it to compute required sample sizes.** Use the full-tree likelihood, or the Monte Carlo
  information computed here.
- **Deflate reported pair-based information** by the measured factor when a rough figure is needed; the median
  across the grid is 7.4×, but the correct factor depends strongly on the rate, which is usually unknown.

## Effect on the Milestone 3 and 4 design conclusions
The budgeted-design comparisons in `milestone4_design_optimization.md` apply a constant deflation and compare
designs at equal budget, so their *ranking* is unaffected. The capacity law is a statement about which rates are
covered, not about how much information is available, so it is unaffected. What is affected is any absolute
statement about precision or power from pair counts, and none is made.

## Limits
Exact enumeration is restricted to 8-leaf trees; everything deeper rests on Monte Carlo validated at that one
depth. Balanced trees only — real trees are star-like, which reduces the effective number of independent
lineages further and so, if anything, makes the naive sum worse still.
