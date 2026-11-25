# Methodology

## The model

Cell state is a binary indicator evolving along a lineage as a two-state continuous-time
Markov chain. Clone `c` has equilibrium occupancy `pi_c` and total switching rate `k_c`.
For two cells of clone `c` whose most recent common ancestor sits at distance `L` from
each,

    Cov(X_i, X_j | c, L) = pi_c (1 - pi_c) exp(-2 k_c L)

(Felsenstein 1985). Averaging over clones and applying the law of total covariance:

    Cov(X_i, X_j | L) = E_c[ pi_c (1-pi_c) exp(-2 k_c L) ]  +  Var_c(pi_c)
                        \________ A(L), state memory ______/   \__ B = V_pi __/

## The two quantities, and why they differ

**Proposition (plateau invariance).** `A(L) -> 0` as `L -> inf` for every distribution of
rates with `k_c > 0`, so the plateau equals `Var_c(pi_c)` exactly. The plateau is a
property of **occupancy alone** and does not move with the switching rate.

This corrects a reading that is easy to fall into: a high plateau is *not* evidence of slow
switching. Two clone populations with identical occupancy distributions and rates differing
by two orders of magnitude have the same plateau.

**The identifiability window.** Writing `u = kL`, the Fisher information about `log k` from
a single pair is proportional to `u^2 exp(-2u) / (1 - exp(-2u))` times an occupancy factor.
On an exact grid it peaks at `u* = 0.43` (varying 0.39–0.49 with occupancy) and keeps at
least half its maximum over `u` in [0.10, 1.0]. A dataset spanning ancestral distances
`[L_min, L_max]` therefore resolves rates only over roughly

    1 / (2 L_max)   <~   k   <~   1 / (2 L_min)

Occupancy carries no such restriction: every distance contributes to the plateau.

**Rate–time confounding.** The likelihood satisfies `L(pi, k; L) = L(pi, k/c; cL)`. With
branch lengths known only up to a global scale, `k` is identified only in units of that
scale. This is why unit edge lengths yield a rate in unknown units, and why a single global
rescaling cannot repair it — the rescaling moves *along* the degeneracy, not across it.
Empirically a naive global scale makes calibration worse, not better.

## Estimation

Per-clone maximum likelihood by Felsenstein pruning over the tree with its given branch
lengths (`src/lineage_inference/ctmc.py`: `VarTree`, `fit_var`; `perclone.py`: `fit_clone`).
Profile-likelihood intervals on a grid in `log k` and `logit pi`.

Variance components subtract estimation noise:

    V_pi    = Var_c(pi_hat_c)    - mean(se^2)
    V_log k = Var_c(log k_hat_c) - mean(se^2)

Raw and noise terms are reported separately in the deposited tables so the correction is
auditable. Shared-`pi` and shared-`k` hypotheses are tested by likelihood ratio against the
per-clone fit.

## Null calibration

The chi-squared reference for the occupancy test is **anti-conservative** on simulated
trees: type-I error 0.143 (95% CI 0.108–0.188, 300 replicates) against a nominal 0.05.
On real trees it errs in both directions, with dataset-specific 95% thresholds from 10.9
to 46.4 against the chi-squared value of 19.68.

All real-data tests therefore use a **parametric bootstrap on the dataset's own trees**:
all clones take the pooled occupancy, each keeps its fitted rate, states are simulated
along the observed topology, and the whole fitting procedure is re-run. 120 draws per
dataset, giving about 0.008 resolution on the p-value.

## The observation chain

Real experiments deliver neither the true genealogy nor true branch times. The simulator
(`src/branching/`) makes each step explicit:

    T0  true genealogy        every cell that ever divided
    T1  sampled genealogy     ancestry among harvested cells
    T2  recorded barcodes     editing, saturation, dropout
    T3  reconstructed tree    greedy / neighbour-joining / UPGMA

By the data-processing inequality, information about `(pi_c, k_c)` is non-increasing along
this chain. Which parameter survives is an empirical question, answered in the experiments.

Two findings matter most. Close relatives — the pairs that carry kinetic information — are
already scarce in the **true sampled genealogy**, a median 1.7% of pairs, before any barcode
is read; reconstruction inherits that geometry rather than creating it. And summing Fisher
information over all pairs as if independent overstates the information in a tree by up to
221-fold, because pairs share ancestry.

## Recorder capacity

An idealised law gives `W = 1 + log10(D_eff)` decades of resolvable rate for `D_eff`
effective divisions. Evaluated against the realised geometry of the six published trees it
**over-predicts by a median of 0.70 decades**: real trees deliver 1.11–1.37 decades where
the law promises 1.54–2.30. The shortfall is attributable to star-likeness.

## Benchmarking against existing methods

Four methods on identical inputs, compared only within admissible estimand groups
(`benchmarks/protocols/estimand_matrix.md`, fixed before any run):

| Method | Native estimand | Branch lengths |
|---|---|---|
| this work | per-clone `(pi_c, k_c)`, variance components, intervals | yes |
| PATH | phylogenetic autocorrelation; rate matrix `Q` | via mean patristic distance |
| PASTRI | per-cell-cycle transition matrix `T`, per depth window | **no** — topological depth |
| SCOUT (OUwie engine) | OU `alpha`, `theta`, `sigma^2` on a continuous trait | yes, unit edges if absent |

PASTRI returns probabilities, not rates: for a two-state chain the off-diagonal sum over
one generation is `1 - exp(-k)`, so it is converted with `k = -log(1 - p)`.

Results are in `benchmarks/reports/benchmark_results.md`. In outline: all four recover
occupancy (correlation 0.81–0.94 with truth; median rank agreement 0.85 on real data),
while the rate estimators disagree by more than an order of magnitude and bracket the truth
from both sides — PATH low by 1.6–2.3×, PASTRI high by 11–16×, this work high by 1.3–1.6×.
Discarding branch times shifts every rate by about four-fold and leaves occupancy
bit-identical.

## References

- Felsenstein J. (1985) Phylogenies and the comparative method. *Am Nat* 125:1–15.
- Simeonov K.P. *et al.* (2021) *Cancer Cell* 39:1150–1162.
- Quinn J.J. *et al.* (2021) *Science* 371:eabc1944.
- Schiffman J.S. *et al.* (2024) *Nat Genet* 56:2174–2184. (PATH)
