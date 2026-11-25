# Estimand compatibility matrix

Written **before** any comparative experiment was run, and not revised afterwards. It
decides which comparisons are admissible. Methods that estimate different quantities are
not placed on a common leaderboard.

Software versions and commits are in `results/tables/TableB2_software.csv`.

---

## 1. What each method natively estimates

### PATH (Schiffman et al.; landau-lab/PATH v1.0, commit `caffc51`)

| Attribute | Value |
|---|---|
| Native estimand | (a) `xcor()`: phylogenetic auto/cross-correlation of states — a Moran's-I style statistic on a phylogenetic weight matrix, with an analytical z-score and p-value. (b) `PATH_inf()`: a cell-state transition probability matrix `P` at t = 1 and rate matrix `Q = logm(P)`. |
| Input phenotype | Categorical (one-hot via `catMat`); continuous also accepted by `xcor`. |
| Tree requirements | Topology. `one_node_depth()` builds the weight matrix from **one-node-depth pairs** (sister cells). Branch lengths enter only through the mean patristic distance `mspd`; `impute_branches=TRUE` instead derives `mspd` from a birth–death–sampling model. |
| Dynamics assumed | Time-homogeneous Markov state transitions; stationarity for the correlation statistic. |
| Output | Point estimates. `xcor` has an analytical z-test. `PATH_inf` has **no** uncertainty interval. |
| Compatible simulations | All of S1–S8 (binary states). |
| Valid Study D comparator | `Q` off-diagonals ↔ our `k`; stationary distribution of `Q` ↔ our `π`; per-clone fits give `V_π`, `V_log k`. `xcor` ↔ our covariance-vs-distance curve. |
| Comparisons that must NOT be made | PATH's `Q` is in units of the mean patristic distance it was handed, so comparing its absolute rate against ours requires both to use the same branch-length convention. No interval-coverage comparison for `PATH_inf` — it reports none. |

### PASTRI (Yang, Li, Yu et al. 2026; yuxiaochen11/PASTRI v0.0.1, commit `67aefd9`)

| Attribute | Value |
|---|---|
| Native estimand | An N×N transition matrix `T`; `T[i,j]` is the probability **per cell cycle** of moving from state i to j. Rows sum to 1. Estimated per depth window `d ≤ fi_depth`. |
| Input phenotype | Categorical cell-type labels. |
| Tree requirements | Topology **only**. Depth is `lca_normalized_height`, computed from node counts (distance-to-root and distance-to-farthest-leaf). **Branch lengths are never read.** |
| Dynamics assumed | Rates may vary with depth — this is the method's stated contribution. Optional hard constraints (`Bound_Matrix`) encode irreversibility. |
| Output | Point estimate per depth window. Uncertainty in the paper comes from a label-shuffling permutation test and gene bootstrap, **not** per-rate intervals. |
| Compatible simulations | All of S1–S8. S5 (environmental induction) is the scenario its stage-specific design targets. |
| Valid Study D comparator | Off-diagonals of `T` ↔ our per-generation switching probability; stationary distribution of `T` ↔ our `π`. Per-clone runs give `V_π`, `V_log k`. |
| Comparisons that must NOT be made | No interval coverage. Because depth is topological, any benchmark that rescales branch times cannot move PASTRI's estimate — that is a prediction to test, not a deficiency to score. |

### SCOUT (Stuart et al. 2025; hrstuart/SCOUT, commit `c72befe`; engine OUwie 3.0.3)

| Attribute | Value |
|---|---|
| Native estimand | Ornstein–Uhlenbeck parameters for a **continuous** trait: `alpha` (rate of reversion to the optimum, i.e. the decay rate of phylogenetic autocorrelation), `sigma^2` (drift), `theta` (regime optimum). Model selection BM vs OU by AICc. |
| Input phenotype | **Continuous** (normalised expression). Not binary. |
| Tree requirements | Topology; branch lengths used if present. SCOUT's own code substitutes unit edges when absent (`R/dropout_validation.R:470-472`, "No edge lengths, replacing with 1s"). |
| Dynamics assumed | OU process; regimes may be assigned per clade. |
| Output | Point estimates plus AICc; model-selection calls. |
| Compatible simulations | The OU trait generated alongside the binary state in every scenario. |
| Valid Study D comparator | **`alpha` ↔ our `k` only as a decay rate of phylogenetic autocorrelation**, since both parameterise `exp(-rate × 2L)`. **`Var_c(theta_c)` ↔ our `V_π`**, since both measure across-clone variance in the trait's central tendency. Held-out trait prediction is a shared observable. |
| Comparisons that must NOT be made | `alpha` is **not** a binary switching rate and must never be plotted on the same axis as `k` as though it were. Its stationary variance `sigma^2/(2 alpha)` has no binary analogue. |

### Study D (this work; Milestone 6 verified version)

| Attribute | Value |
|---|---|
| Native estimand | Per clone `(π_c, k_c)` by maximum likelihood (Felsenstein pruning); across-clone variance components `V_π` and `V_log k`; likelihood-ratio tests for shared π and shared k; profile-likelihood intervals; an identifiability window. |
| Input phenotype | Binary. |
| Tree requirements | Topology and per-edge lengths (`VarTree`), or unit edges (`FastTree`). |
| Dynamics assumed | Two-state time-homogeneous CTMC per clone. |
| Output | Point estimates, profile intervals, calibrated tests, and a statement of which rates the tree can resolve. |

---

## 2. Admissible comparison groups

**Group A — direct estimand comparison (same quantity, same units).**
Study D `k` vs PATH `Q` off-diagonal vs PASTRI `T` off-diagonal, **all on the same tree with
the same branch-length convention**. All three are per-generation state-change rates for a
binary trait. SCOUT is excluded.

**Group B — occupancy / central tendency.**
Study D `V_π` vs PATH stationary distribution of `Q` vs PASTRI stationary distribution of
`T` vs SCOUT `Var_c(theta_c)`. All four measure across-clone variation in where the trait
sits. SCOUT enters here on its continuous trait, the others on the binary state.

**Group C — evolutionary memory / phylogenetic autocorrelation.**
Study D's covariance-vs-distance curve vs PATH `xcor` vs SCOUT `alpha`. A shared decay-rate
estimand, not a shared parameter value.

**Group D — shared observable prediction.**
Held-out leaf state (binary: Study D, PATH, PASTRI) and held-out leaf trait (continuous:
SCOUT). Measured as mean absolute error on cells withheld from fitting. This is the only
comparison in which all four methods appear on one axis, and it is legitimate because the
quantity scored is a prediction about data, not a parameter.

**Group E — assumption sensitivity.**
How much each method's own estimate moves when the tree geometry, branch-time convention or
reconstruction changes. Each method is compared **against itself**, so differing estimands
do not matter.

**Group F — inferential calibration.**
Empirical coverage of nominal intervals. **Only Study D and SCOUT qualify** (profile
likelihood; OUwie's likelihood surface). PATH `PATH_inf` and PASTRI report no intervals and
are scored as "not applicable", never as zero coverage.

---

## 3. Comparisons explicitly ruled out

1. A single accuracy leaderboard across all four methods.
2. SCOUT `alpha` against Study D / PATH / PASTRI rates on a shared numeric axis.
3. Runtime as a quality measure.
4. Assigning PATH or PASTRI a coverage of zero for not reporting intervals.
5. Interpreting PASTRI's invariance to branch-length rescaling as an error — it is a
   consequence of a topological depth definition and is reported as such.
6. Treating agreement between methods applied to the same dataset as independent biological
   replication.
