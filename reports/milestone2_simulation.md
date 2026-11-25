# Milestone 2: simulation results

Simulator: `src/branching/simulator.py` — an individual-cell birth–death–switching–dissemination process with
clone-specific switching propensity, organ-specific induction, state-dependent dissemination, stochastic founding,
evolving-barcode edits (so only some branches are visible), partial sampling, state misclassification and clone
label loss. Ground truth retained: full ancestry, founder identities and states, switch counts, clone parameters.

Experiments: E8 (`run_e8.py`, mechanism signatures) and E9 (`run_e9.py`, recovery and stress), both executed on the
server. Seeds 20261124 / 20261125, 15–20 replicates per mechanism.

## Mechanism signatures (E8/E9, primary tumour, 15 replicates each)

| Mechanism (truth) | true B | naive statistic: variance of clone fractions | decomposition B̂ | P(CI excludes 0) |
|---|---|---|---|---|
| A stable-state selection | 0 | **0.040** | 0.000 | **0.00** |
| B heritable propensity | 0.078 | 0.068 | 0.061 | 0.00 |
| C organ induction | 0 | 0.002 | 0.000 | **0.00** |
| D founder bottleneck | 0 | 0.002 | 0.000 | **0.00** |
| E coupled (B + C) | 0.067 | 0.066 | 0.055 | 0.40 |

**The headline.** Across 45 replicates of the three mechanisms whose true B is zero, the decomposition produced
**no false positives**. The naive statistic — the variance of clone-level state fractions, which is what a clone
analysis reports by default — reads 0.040 under pure stable-state selection, which is two thirds of its reading
under genuine heritable propensity (0.068). The two are not distinguishable without the decomposition.

**The cost is power.** Under true heritable propensity the point estimate is close to the truth (0.061 against
0.078) but the clone bootstrap CI excluded zero in 0 of 15 replicates, and in 6 of 15 for the coupled mechanism.
With ~20 clones and a few thousand cells, this estimator is **specific but not sensitive**. A null result is weak
evidence; a positive result is strong.

## Observed covariance curves (E9, `results/figures/m2_fig8_mechanism_signatures.png`)
- Clone-centred curves (which remove the plateau by construction) rise with MRCA depth in every mechanism — this
  is state memory, and it is what identifies the decay rate.
- Site-centred curves start near zero at depth 0 for induction (−0.0004) and bottleneck (0.0002), but at 0.058 for
  heritable propensity. That non-zero intercept at maximum lineage distance *is* the plateau.
- Stable-state selection starts at 0.0027 and rises steeply to 0.12: all memory, no plateau.

## Metastatic compartments are uninformative for B
| Mechanism | true B | B̂ in metastases | median clones per site |
|---|---|---|---|
| B heritable propensity | 0.078 | 0.012 | 7 |
| E coupled | 0.067 | 0.023 | 7 |

Metastases are founded by few cells and carry few clones, so B̂ is attenuated by a factor of 3–6. This is
Proposition T4 in practice, and it is why the real-data analysis estimates B in the primary tumour.

## Stress tests (E9, mechanism B, true B = 0.065)
| State misclassification | Clone labels lost | Cells sampled | B̂ |
|---|---|---|---|
| 0 | 0 | 4000 | 0.045 |
| 0.05 | 0 | 4000 | 0.046 |
| 0.15 | 0 | 4000 | 0.019 |
| 0 | 30% | 4000 | 0.017 |
| 0.15 | 30% | 800 | 0.011 |

Both nuisances **attenuate** B̂ and never inflate it. Measurement error therefore biases this analysis toward
concluding "no heritable propensity", which is the safe direction and reinforces that a null is weak evidence.

## Model recovery is not evidence of mechanism (counterexample C8)
A nested-model comparison that selects clone-specific stationary fractions (M1) over a single-state model (M0)
would also do so under pure stable-state selection, because clones genuinely differ in composition in both cases.
The decomposition is the stricter test. Both are reported for the real data, and they disagree — see
`simeonov_reanalysis.md`.
