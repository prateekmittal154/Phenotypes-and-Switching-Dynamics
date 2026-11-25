# Milestone 6, Experiments 1 and 2: reconstruction-aware interval calibration and branch-time models

E20 (`experiments/reconstruction/run_e20.py`), 480 replicates: 2 recorders × 4 true rates × 60 replicates,
96 sampled cells each. Profile-likelihood intervals for the switching rate at six inference variants on
identical observed states. Theory: Proposition V4 in `theory/information_hierarchy.md`.

| Variant | What it conditions on |
|---|---|
| V1 | true sampled topology with **true branch times** (reference) |
| V2 | true topology, unit edges |
| V3 | reconstructed topology, unit edges — **standard practice** |
| V4a | reconstructed, all edges scaled by a single factor $c = D/h$ |
| V4 | reconstructed, **ultrametric** branch times (every leaf at age $D$) |
| V5 | ensemble over three reconstructions, ultrametric |
| V6 | V5 integrating over an **uncertain** time scale |

## Result 1 — the inference machinery is correct when ancestry and time are known
V1 attains **0.927 [0.900, 0.947]** coverage against a nominal 0.95, with a median log bias of +0.124. Every
subsequent loss is therefore attributable to the observation process, not to the estimator.

## Result 2 — pooled calibration across all regimes

| Variant | Coverage [95% Wilson] | median log bias | median log interval width |
|---|---|---|---|
| V1 | **0.927 [0.900, 0.947]** | +0.124 | 4.13 |
| V2 | 0.542 [0.497, 0.586] | +1.070 | 3.61 |
| V3 | 0.494 [0.449, 0.538] | +1.328 | 3.61 |
| **V4a** | **0.362 [0.321, 0.406]** | **+1.975** | 3.10 |
| V4 | 0.504 [0.460, 0.549] | +1.413 | 3.61 |
| V5 | 0.600 [0.556, 0.643] | +1.070 | 3.87 |
| **V6** | **0.675 [0.632, 0.715]** | +1.314 | 4.13 |

**A failed correction, reported because it is informative.** Rescaling every edge by a single factor $D/h$ —
the obvious first attempt — makes calibration *worse* than doing nothing (0.362 against 0.494) and inflates the
bias to +1.975. Proposition V4 explains why: a global scale makes only the deepest root-to-leaf path equal the
elapsed time, so in an unbalanced tree the many shallower leaves get paths that are too short.

**Reconstruction-aware inference does improve calibration**: 0.494 → 0.600 (ensemble) → **0.675** (ensemble
with an uncertain time scale), a 37% relative gain, achieved by widening intervals to reflect the uncertainty
that was previously ignored rather than by sharpening estimates.

## Result 3 — the gain is concentrated exactly where information exists

Pooling across regimes understates the result badly. Coverage by true rate:

| true k | regime | V1 | V3 (standard) | V5 | V6 |
|---|---|---|---|---|---|
| 0.05 | **below** the window | 0.842 | 0.083 | 0.292 | 0.342 |
| 0.30 | below / at the edge | 0.892 | 0.033 | 0.167 | 0.375 |
| **1.00** | **within** the window | 0.983 | 0.883 | **0.967** | **1.000** |
| 3.00 | above the window | 0.992 | 0.975 | 0.975 | 0.983 |

**Within the identifiable window, reconstruction-aware inference is essentially correctly calibrated** (V5
0.967, V6 1.000, against V3's 0.883). Below the window every variant fails, V1 excepted — coverage of 3–38%
even with reconstruction uncertainty propagated.

This is the precise answer to the question the brief posed. Reconstruction-aware inference improves
calibration where the recorder wrote the information, and **does not manufacture coverage where it did not**.
That is the behaviour Proposition V1 requires: no processing step can create information, so a method that
reported nominal coverage below the window would be wrong, not impressive.

A caveat on the k = 3 row: coverage near 1.00 there comes with log interval widths of 3.6–4.1, i.e. intervals
spanning a 35- to 60-fold range. Coverage without precision is not usable inference, and these intervals say
little more than "the rate is fast".

## Result 4 — branch-time calibration (Experiment 2)
Of the 1.83-fold bias Milestone 5 attributed to the clock assumption:

- **A single global scale removes none of it and adds more** (+1.33 → +1.98 in log units).
- **Ultrametric calibration**, which uses only the fact that all sampled cells are the same age, is the valid
  correction in principle but reduces bias only modestly here (+1.33 → +1.41 on the pooled median; it helps at
  k = 0.05, +2.22 → +1.89).
- **Ensembling plus ultrametric times** gives the best bias (+1.07, matching V2's true-topology bias).
- **Integrating over an uncertain time scale** buys coverage rather than bias: +1.31, but coverage 0.675.

**Edit distance is not elapsed time, and we did not impose an external clock.** The ultrametric construction
uses only the experiment's own knowledge that all harvested cells are the same age; no edit rate is assumed.
The residual bias is irreducible without an independent clock, which is why the design report recommends
measuring the edit rate rather than inferring time from edits.

## Limits
One tree size (96 cells), one generation count, synchronous divisions. Profile intervals assume the usual
$\chi^2_1$ cut-off for a single parameter, which V1's 0.927 coverage validates in the correctly specified case.
The uncertain-scale prior (log-normal, sd 0.3) is a stand-in, not elicited from data.
