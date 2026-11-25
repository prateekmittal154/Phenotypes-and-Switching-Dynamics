# Benchmark design

Fixed before results were inspected. Deviations are logged in `reports/deviations.md`.

## Resource budget, set in advance

A single analysis server (24 cores, shared with other users). Each method gets one core per
job, with BLAS threading disabled so the job count equals the core count. No method is given
more compute than another on the same replicate. Replicate counts were chosen so every
benchmark finishes inside the budget, and were not revised after seeing any comparison.

## Benchmark I — ground truth, eight mechanisms

200 replicates: 8 scenarios x 25. Each replicate has 12 clones x 48 cells = 576 cells,
11 generations of synchronous division.

| Scenario | Generating mechanism | True `V_pi` | True `V_log k` |
|---|---|---|---|
| S1 | stable-state selection: occupancy differs, k = 0.02 | > 0 | 0 |
| S2 | equal occupancy (pi = 0.5), k varies over 0.05-2.0 | 0 | > 0 |
| S3 | occupancy varies, shared k = 0.3 | > 0 | 0 |
| S4 | both vary | > 0 | > 0 |
| S5 | environmental induction: k rises 4-fold after generation 5 | > 0 | 0 |
| S6 | founder bottleneck at generation 4, 6 survivors | > 0 | 0 |
| S7 | unequal clonal expansion, per-clone division probability 0.55-1.0 | > 0 | 0 |
| S8 | observation misspecification: 10% state flips, trait noise SD 0.5 | > 0 | 0 |

S2 and S3 are the discriminating pair. A method that recovers occupancy but not kinetics
should report structure in S3 and little in S2; a method that confuses the two will report
heterogeneity in both.

The continuous trait is an Ornstein-Uhlenbeck process on the same genealogy, with clone
optimum `theta_c = 2(pi_c - 0.5)` and reversion rate `alpha_c = k_c`. This gives SCOUT a
native-format input instead of a binary variable it was not designed for, while keeping the
generating structure aligned so the comparison means something.

## Benchmark II — lineage distance and sampling

75 replicates: 5 geometries x 15. The true rate is **k = 0.3 for every clone in every
geometry**, so the true `V_log k` is exactly zero and any reported rate heterogeneity is an
artefact of geometry.

| Geometry | Construction |
|---|---|
| G_shallow | 7 generations, random sampling |
| G_deep | 13 generations, random sampling |
| G_clustered | 11 generations, cells drawn from whole clades |
| G_dispersed | 11 generations, cells spread across distinct clades |
| G_starlike | 11 generations with a bottleneck at generation 3 |

For each replicate we record the realised within-clone ancestral-distance distribution and
the fraction of pairs whose `kL` falls in the informative band [0.10, 1.0]. That fraction is
the independent variable: it tests whether the theoretical window predicts actual
performance.

## Benchmark III — branch-time calibration

Every Benchmark I replicate is analysed twice: once with true branch lengths in generations,
once with every edge set to 1. Nothing else changes, so the difference isolates the clock
assumption from everything else. Each method is compared against itself, so differing
estimands do not matter.

## Benchmark IV — reconstruction uncertainty

20 replicates, 8 clones x 48 cells. Barcodes are simulated on the full genealogy with a
Cassiopeia-like recorder (30 sites, per-site edit rate 0.06 per generation, 40 alleles, 10%
dropout, 10% missing) and trees are rebuilt by greedy, neighbour joining and UPGMA. Five
trees per replicate over identical cells and states: true-with-times, true-with-unit-edges,
and the three reconstructions. Methods that take a single tree are run separately on each
and the spread reported; that is **not** called Bayesian uncertainty integration.

## Benchmark V — real data

116 clone trees from two published systems, exported through the same loaders the canonical
Study D analysis uses, so the trees and state definitions are identical.

| Animal | Clone trees | Cells |
|---|---|---|
| Simeonov M1 | 8 | 615 |
| Simeonov M2 | 21 | 1,706 |
| Quinn 5k | 64 | 5,999 |
| Quinn 10k | 9 | 801 |
| Quinn 30k | 1 | 132 |
| Quinn 100k | 13 | 1,141 |

Animals are analysed separately and never pooled. The Quinn 30k animal yields a single
eligible tree, so no across-clone variance can be formed from it; it is carried through the
tables as a single clone and excluded from any across-clone statistic.

Continuous traits: the published EMT programme score for Simeonov, the cell-cycle score for
Quinn. Binary states: the same thresholded definitions the canonical analysis uses.

## Uncertainty reporting

- Coverage and other proportions: Wilson intervals.
- Medians and biases: percentile bootstrap, 2,000 resamples.
- Interval coverage is evaluated **only** for methods that report intervals (Study D,
  SCOUT). PATH's `PATH_inf` and PASTRI report none and are recorded as "not applicable".

## Settings chosen for the competing methods

- **PATH**: `PATH_inf(tree, nstates = 2)` with the package default weight matrix
  (`one_node_depth`) and `impute_branches = FALSE`, so the mean patristic distance comes
  from the tree supplied. `xcor` is used with the same weight matrix.
- **PASTRI**: `fi_depth` is swept over {0.2, 0.4, 0.6, 0.8, 1.0} and the headline value is
  0.4, matching the cumulative window the authors use for their C. elegans benchmark. A
  single fixed 1.0 would include the distal adjacencies the authors themselves report cause
  overestimation. `Bound_Matrix` imposes no irreversibility, since the two-state process is
  reversible by construction.
- **SCOUT**: the OUwie engine with `model = "OU1"` and `"BM1"`, `scaleHeight = FALSE`,
  one regime per clone. Unit edges are substituted when a tree carries none, which is
  SCOUT's own rule (`R/dropout_validation.R:470-472`).
- **Study D**: `VarTree` + `fit_var`, the Milestone 6 verified estimator, with per-edge
  lengths when supplied.

Conversion to a common scale: PASTRI returns a per-cell-cycle transition **probability**
matrix, so its off-diagonal sum `p` is converted to a rate by `k = -log(1 - p)`. PATH's `Q`
is already a rate. Study D's `k` is the total rate of the two-state generator. All three are
then per generation.
