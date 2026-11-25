# Claim-to-evidence audit

Every claim the manuscript makes, the evidence behind it, and its class. Claims that fail the audit are
removed from the manuscript rather than softened.

**Evidence classes:** *identity* (algebra, given stated assumptions) · *numerical* (deterministic computation)
· *simulation* (known ground truth, with Monte Carlo uncertainty) · *observational* (real data, no ground
truth).

## Mathematical claims

| Claim | Evidence | Class | Verdict |
|---|---|---|---|
| Within-clone pair covariance = state memory + occupancy plateau | Lemma T1 + law of total covariance; proofs in `theory/proofs/T1_T2.md` | identity | **holds**; Lemma T1 is classical (Felsenstein 1985) and is attributed as such |
| The plateau is invariant to the switching rate | Proposition U1 | identity | **holds**; corrects the Milestone 2 reading |
| Induction rescales but cannot create a plateau | Proposition U3 | identity | **holds** in the stated weak form only |
| $(\pi_c, k_c)$ identified iff $1/(2L_{\max}) \lesssim k \lesssim 1/(2L_{\min})$ | Proposition U4 | identity | **holds** |
| $kL^\ast \approx 0.43$, flat over a ten-fold band | E12, exact grid | numerical | **holds**; replaces the Milestone 3 value 0.405 |
| Capacity law $W = 1 + \log_{10} D_{\text{eff}}$ | Steps 1,3,4 identities on a numerical band | identity + numerical | **holds as an idealisation**; must be presented with E16 |
| Information-processing chain across T0–T3 | Proposition V1 | identity | **holds** |
| Rate–time confounding: $\mathcal{L}(\pi,k;L) = \mathcal{L}(\pi,k/c;cL)$ | Proposition V4 | identity | **holds**; explains the clock bias |

## Numerical and simulation claims

| Claim | Evidence | Verdict |
|---|---|---|
| Naive pair-summing overstates information 1–221× | E13, exact enumeration (8 leaves) + MC validated to 6.0% | **holds** |
| χ² reference for the occupancy test is anti-conservative (0.143 [0.108, 0.188]) | E14, 300 replicates, Wilson | **holds** |
| Occupancy-variance interval coverage collapses to 0.253 below the window | E14 | **holds** |
| Star-likeness is present in the *true* sampled genealogy (1.6–1.9% close pairs) | E17, 480 replicates | **holds** — the central Milestone 5 claim |
| Clock miscalibration costs 1.83×; topology error adds 1.04–1.10× | E17, paired within-replicate | **holds** |
| Recorder capacity doubles topology recovery | E17 | **holds**; greedy result is *not* a verdict on Cassiopeia |
| Inference is correctly calibrated given true ancestry and times (0.927) | E20, 480 replicates | **holds** |
| A single global scale correction makes calibration worse (0.494 → 0.362) | E20 | **holds** — a documented failed correction |
| Reconstruction-aware inference raises coverage 0.494 → 0.675 pooled | E20 | **holds** |
| Within the window it reaches 0.967–1.000; below it, 3–38% | E20 by regime | **holds** — the key qualification |
| Two-timepoint sampling triples informative pairs | E21, 40 observations/strategy | **holds** for the pair fraction |
| Two-timepoint sampling improves coverage 0.725 → 0.825 | E21 | **suggestive only** — Wilson intervals overlap; labelled as such |

## Observational claims (real data)

| Claim | Evidence | Verdict |
|---|---|---|
| Clone-associated occupancy is heterogeneous | E22, dataset-specific bootstrap, 8 of 8 eligible datasets, 5 animals | **holds** |
| Quinn 5k primary-only is borderline (p = 0.050) | E22 | **holds**, and is reported as borderline |
| Occupancy is stable across reconstructions (±6% Quinn, 3 d.p. Simeonov) | E18, E19 | **holds** |
| Rate heterogeneity varies 6.5-fold across published reconstructions | E18 | **holds** |
| Metastatic EMT enrichment is compositional | M2 decomposition; composition term uses primary expression only | **holds** |
| Within-clone cross-tissue shifts are depth-confounded | M2 and M4 controls | **holds**; the within-clone channel is explicitly not interpreted |
| Switching rates are unidentifiable in both systems | M4 window + E18 spread + E20 coverage | **holds** |

## Claims that were removed
| Removed claim | Why |
|---|---|
| The plateau measures "heritable switching propensity" | Invariant to the rate (U1) |
| $kL^\ast = 0.405$ is constant | Varies 0.39–0.49 with occupancy |
| Recorders span 0.09–0.24 decades | Arbitrary criterion; 1.1–1.4 on an information basis |
| Cassiopeia depth raises the usable ceiling to 0.86 | False on the realised geometry; its empirical ceilings are lower |
| "0% false positives, 100% power" | From 10 replicates; the occupancy test's true type-I error is 14% |
| Six eligible animals | Under Milestone 6's stricter tree caps the correct count is 5 |

## Claims deliberately not made
Any mechanistic statement about why clones differ in occupancy; any human relevance; any claim that better
reconstruction could recover the rates; any gene-level or spatial validation; any cost-per-precision figure
(the pair-based power calculation is unreliable by 1–221×).
