# Milestone 3: lineage sampling design

Milestone 2 noticed that the spread of observable MRCA depths seemed to matter more than the number of cells
sequenced. This report formalises that, and then tests it against two real lineage-recording technologies.

## 1. The information calculus (Proposition U5)

For a within-clone binary pair at lineage distance L, with occupancy π and switching rate k, write
ρ = e^{−2kL}. The pair falls into one of four cells with probabilities p₁₁ = π² + πqρ, p₀₀ = q² + πqρ,
p₁₀ = p₀₁ = πq(1 − ρ). The Fisher information about log k is

  I_log k(L) = I_ρ(ρ) · (2 k L ρ)²,  I_ρ(ρ) = (πq)² [1/p₁₁ + 1/p₀₀ + 2/p₀₁].

**The optimum is interior and essentially universal.** Information vanishes as L → 0 (relatives too close, ρ ≈ 1,
no contrast) and as L → ∞ (relatives too distant, ρ ≈ 0, nothing left to measure). Numerically the maximum sits at

| k | π | optimal L | k·L\* |
|---|---|---|---|
| 0.1 | 0.3 | 4.07 | 0.407 |
| 0.1 | 0.5 | 4.07 | 0.407 |
| 0.5 | 0.3 | 0.81 | 0.405 |
| 0.5 | 0.5 | 0.81 | 0.405 |
| 2.0 | 0.3 | 0.21 | 0.421 |
| 2.0 | 0.5 | 0.20 | 0.394 |

**k L\* ≈ 0.405**, stable across two decades of k and across occupancy. The simplified maximiser of L²e^{−4kL}
gives k L\* = 0.5, which is the right order and the easier rule of thumb. In words: **the most informative pairs
are those separated by roughly one switching time.**

Information about occupancy π behaves completely differently: it comes from the leaf marginal and grows with the
**number of cells**, nearly independently of depth.

### The design rule
> **Occupancy is bought with cells. Switching rate is bought with pairs separated by about one switching time.**

Sequencing more cells whose pairwise lineage distances all sit far from 0.4/k buys almost nothing about dynamics.
This is a falsifiable, quantitative recommendation, and it inverts the usual instinct to scale cell number.

> **Superseded in part by Milestone 4.** The window figures below use an ad-hoc criterion
> (1/(2L_max), 1/(2L_min)). On a Fisher-information criterion (>=50% of peak information, k*L in
> [0.10, 1.0]) the same recorders cover **1.1-1.2 decades**, not 0.09-0.24. The qualitative
> conclusion is unchanged -- the fitted Simeonov rates still sit above the covered range -- but the
> quantitative claim was too pessimistic. See `reports/milestone4_design_optimization.md`.

## 2. Do real recorders satisfy the window?

Proposition U4 gives the identifiable range of k as 1/(2 L_max) ≲ k ≲ 1/(2 L_min), set by the spread of observed
MRCA depths. Computed directly from the published trees of two different technologies
(`results/tables/m3_depth_spectra_summary.csv`, median over clones):

| System | Dataset | leaves/clone | max depth | k identifiable from | to | window (decades) |
|---|---|---|---|---|---|---|
| macsGESTALT | mouse 1 | 74 | 4.0 | 0.167 | 0.250 | **0.18** |
| macsGESTALT | mouse 2 | 64 | 3.5 | 0.250 | 0.500 | **0.18** |
| Cassiopeia | m5k | 76 | 14.0 | 0.038 | 0.071 | **0.24** |
| Cassiopeia | m10k | 93 | 8.0 | 0.071 | 0.100 | **0.18** |
| Cassiopeia | m30k | 400 | 13.0 | 0.042 | 0.051 | **0.09** |
| Cassiopeia | m100k | 35 | 9.0 | 0.062 | 0.125 | **0.20** |

**Neither technology spans even a quarter of one decade.** A rate can be identified only if it already happens to
lie inside a window about 1.3–1.7-fold wide. Cassiopeia's deeper trees (max depth 8–14 against 3.5–4) shift the
window to *slower* rates, but do not widen it.

This is the single most generalisable finding of Milestone 3: **the identifiability limitation found in
Milestone 2 is not a quirk of macsGESTALT. It is a shared property of current lineage recorders**, and no amount
of extra sequencing fixes it, because the binding constraint is the depth spectrum rather than the cell count.

## 3. Equal-budget designs: the trade-off is real and quantified

Every design below spends the **same 1,280 cells**, split differently between clones and tree depth
(`results/tables/e11_equal_budget.csv`, 10 replicates, medians). In this simulator a balanced tree of depth d has
2^d leaves, so depth and cells-per-clone move together; the comparison is therefore "few deep clones versus many
shallow clones at fixed total cells", which is the choice an experimenter actually faces.

| Design | clones | cells/clone | depth | true V_π | estimated V_π | true V_log k | estimated V_log k | rates off-bound |
|---|---|---|---|---|---|---|---|---|
| many clones, shallow | 40 | 32 | 5 | 0.030 | 0.041 | 0.859 | 0.377 (44%) | 31/40 |
| balanced | 20 | 64 | 6 | 0.032 | 0.037 | 0.633 | 0.070 (11%) | 18/20 |
| **few clones, deep** | 10 | 128 | 7 | 0.034 | **0.034** | 0.905 | **0.747 (83%)** | 9/10 |
| very few, very deep | 5 | 256 | 8 | 0.028 | 0.021 | 0.558 | 0.247 (44%) | 5/5 |

**The two components want opposite designs.**
- **Occupancy V_π** is recovered well by every design with ≥ 10 clones and is *underestimated* with only 5
  (0.021 against 0.028): it needs **clones**, because it is a variance across clones.
- **Switching heterogeneity V_log k** needs **depth**: the 10-clone deep design recovers 83% of it, the
  40-clone shallow design 44%, and the balanced design only 11%.
- Pushing to 5 very deep clones is worse for both, because a variance over 5 units is itself poorly estimated.

At this budget the best compromise is roughly **10–20 clones sampled deeply** rather than the instinctive "as
many cells as possible, spread widely". Note also that these are recovery fractions, not coverage: calibrated
interval coverage for V_log k was not established, and that remains a gap.

## 4. What would fix it
Widening the window requires edits that accumulate over a wide range of lineage depths, not more cells at the
depths already sampled. Concretely:
1. **Stagger the recording.** Induce editing in several waves, or use recorders with heterogeneous edit rates, so
   that shallow and deep coalescences are both resolved. A recorder with a single fast burst saturates and
   produces the narrow spectra above.
2. **Target the expected timescale.** Given a prior guess for k, place the bulk of resolvable depths near
   L = 0.4/k. If k is unknown within an order of magnitude, deliberately spread edit rates across that range
   instead of optimising for one.
3. **Spend the marginal dollar on depth diversity, not cells.** Equal-budget simulations are in
   `results/tables/e11_equal_budget.csv`.
4. **Report the window.** Any paper estimating switching rates from lineage trees should state
   1/(2L_max) and 1/(2L_min) for its data, so readers can see whether the estimate was identifiable at all. None
   of the methods surveyed does this.

## 5. Limits of this calculus
- It is derived for a two-state, time-homogeneous, stationary process on a tree with known topology. Lineage
  reconstruction error, non-stationarity and state misclassification all reduce information further, so the
  windows above are optimistic.
- Branch lengths are taken as the number of intervening edits, which is a clock only if edits accumulate at a
  constant rate. Where that fails, L is a monotone but non-linear function of time, and the windows shift.
- The optimum k L\* ≈ 0.405 assumes the pair covariance is the only statistic used. A full-tree likelihood uses
  more, so the true optimum is somewhat broader; the simulation comparison in E11 addresses this empirically.
