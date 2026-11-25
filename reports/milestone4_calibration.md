# Milestone 4: calibration of tests and intervals

300 replicates per scenario, 12 clones each, parametric-bootstrap null distributions, Wilson intervals on every
rate. Script `experiments/calibration/run_e14.py`; tables `e14_calibration_summary.csv`, `e14_thresholds.csv`.

## Scenarios
| Label | Truth | What it tests |
|---|---|---|
| H0_pi | occupancy shared, rates vary | type-I error of the **shared-occupancy** test |
| H0_k | occupancies vary, rate shared | type-I error of the **shared-rate** test |
| H0_k_frozen | occupancies vary, k = 0.01 (below the window) | the boundary case |
| H1_pi / H1_k | the corresponding alternatives | power |
| MIS | three-state generator collapsed to two, all clones identical | out-of-family misspecification |

## Result 1 — the χ² reference is wrong in both directions, and it matters which

| Test | Scenario where the null holds | Rejection rate (95% Wilson) | Verdict |
|---|---|---|---|
| shared-occupancy, χ² | H0_pi | **0.143 [0.108, 0.188]** | **nearly 3× inflated** |
| shared-occupancy, χ² | MIS (clones identical) | **0.123 [0.091, 0.165]** | inflated under misspecification too |
| shared-rate, χ² | H0_k | **0.013 [0.005, 0.034]** | conservative |
| shared-rate, χ² | H0_k_frozen | **0.010 [0.003, 0.029]** | conservative at the boundary |

The bootstrap thresholds show the same thing directly: for the occupancy test the calibrated 95% cut-off is
**32.59** against the χ² value of 19.68, while for the rate test it is **17.35**, slightly *below* χ².

**This corrects a Milestone 3 claim.** From 10 replicates I reported "0% false positives". That is right for the
rate test — the authoritative figure is 1.0–1.3%, genuinely conservative — but wrong for the occupancy test,
whose true type-I error is about 14%, not 5%.

**Do the real-data occupancy conclusions survive?** Yes, with room to spare. Against the calibrated threshold of
32.59 rather than χ²:

| Dataset | occupancy LRT statistics | exceed calibrated threshold |
|---|---|---|
| Simeonov (2 mice, 2 scopes) | 106.6, 74.9, 93.0, 145.8 | **4 of 4** |
| Quinn, all sites (4 mice) | 164.0, 7.0, 56.0, 162.3 | **3 of 4** (the 3-clone 30k animal fails) |
| Quinn, primary only, depth-adjusted | 45.9, 59.3, 89.0 | **3 of 3** |

The statistics are three- to five-fold above the calibrated cut-off, so a threefold inflation of the reference
does not threaten them.

## Result 2 — power

| Test | Scenario | Power at the calibrated threshold |
|---|---|---|
| shared-rate | H1_k | **0.883 [0.843, 0.914]** |
| shared-occupancy | H1_pi | **0.493 [0.437, 0.549]** |

The rate test is well powered when rates genuinely differ. The occupancy test has only about 50% power at 12
clones once the threshold is calibrated honestly — so a null occupancy result is weak evidence, while a positive
one is strong.

## Result 3 — interval coverage, including the regime where it collapses

Coverage of the 95% clone-bootstrap interval for the occupancy variance $V_\pi$:

| Scenario | Coverage (95% Wilson) |
|---|---|
| H0_pi | 0.967 [0.940, 0.982] |
| H0_k | 0.920 [0.884, 0.946] |
| H1_k | 0.947 [0.915, 0.967] |
| H1_pi | 0.880 [0.838, 0.912] |
| MIS | 0.990 [0.971, 0.997] |
| **H0_k_frozen** | **0.253 [0.207, 0.305]** |

Coverage is acceptable everywhere except the frozen regime, where it collapses to a quarter of nominal. This is
the Milestone 3 failure mode confirmed at a serious sample size: **below the identifiability window the interval
is not merely wide, it is wrong.**

Per-clone profile intervals for $\pi_c$ cover at 0.91–0.97 in every scenario except misspecification, where they
fall to **0.719** — a three-state process read as two states breaks the per-clone intervals even though the
aggregate $V_\pi$ interval survives.

## What to do about it, and what the analysis already does
1. **Do not use χ² for the occupancy test.** Use the parametric-bootstrap threshold. The real-data conclusions
   were re-checked against it above.
2. **The rate test may be used as reported**; it is conservative, which is the safe direction for a claim we are
   not making.
3. **Report the bound flag and refuse point estimates below the window.** Interval coverage there is 25%.
4. **Treat null occupancy results as uninformative** at this clone count, given 49% power.

## Monte Carlo uncertainty
Every rate above is a proportion out of 300 with its Wilson interval quoted. The widest interval is about
±0.05, so distinctions smaller than that are not resolved by this experiment.

## Limits
One clone count (12), one tree size, one misspecification family (a collapsed three-state chain). Branch-length
error, reconstruction error and non-Markov dynamics are not covered here; reconstruction error is addressed
separately in Milestone 5.
