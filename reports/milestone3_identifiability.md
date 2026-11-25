# Milestone 3: identifiability and estimator comparison

## 1. The question, restated correctly
Milestone 2 asked whether a clone's lineage carries "heritable switching propensity". The audit (U1) shows that
question was conflated. Two separable quantities exist:

- **V_π = Var_c(π_c)** — clone-associated *equilibrium occupancy*: how much clones differ in how often they sit
  in the high state. This is what the Milestone 2 plateau measured.
- **V_logk = Var_c(log k_c)** — clone-associated *switching dynamics*: how much clones differ in how fast they
  interconvert. This is what "heritable switching propensity" should mean, and Milestone 2 never estimated it,
  because it fitted a single shared rate.

A clone can have high occupancy because it switches into the state rapidly, or because it was born there and
never leaves. Those are different biology and the plateau cannot tell them apart.

## 2. Necessary and sufficient conditions (U4)
For one clone, observed at stationarity, with leaf pairs at two or more distinct lineage distances:
- π_c is identified from the leaf marginal;
- k_c is identified from the pair covariance, ρ(L) = e^{−2k_cL}, **provided** k_c lies inside
  **1/(2 L_max) ≲ k_c ≲ 1/(2 L_min)**;
- two or more distinct distances over-identify the model and make the two-state assumption testable.

Outside the window the likelihood is flat (k indistinguishable from 0) or saturated (indistinguishable from ∞).
The per-clone estimator reports an explicit `k_at_bound` flag whenever this happens, and such clones are excluded
from V_logk rather than silently contributing a bound value.

## 3. Estimators compared
| # | Estimator | What it targets | Model dependence |
|---|---|---|---|
| 1 | Milestone 2 two-step covariance decomposition | pooled plateau = V_π only | low; assumes a shared k |
| 2 | **Per-clone ML with profile intervals + noise-corrected variance components** (new) | π_c, k_c, then V_π and V_logk | medium; two-state CTMC per clone |
| 3 | Likelihood-ratio tests of a shared π or shared k across clones | significance of each component | medium |
| 4 | Profile-likelihood intervals per clone, with a bound flag | partial identification when the window fails | medium, but honest about failure |
| 5 | PATH-style single phylogenetic autocorrelation | one population-average statistic | low; cannot separate the two components |

Estimator 2 subtracts the mean squared standard error from the raw variance of the per-clone estimates, so that
estimation noise is not read as between-clone heterogeneity — the same method-of-moments correction that
counterexample C6 showed to be necessary.

**Does the full tree beat aggregated pairwise covariance?** Yes for occupancy and, more importantly, the full
tree is what makes per-clone rates estimable at all: the pooled covariance curve averages over clones and so can
only return one rate. The cost is model dependence — the per-clone likelihood commits to a two-state
time-homogeneous CTMC with known topology, whereas the pairwise curve is nearly assumption-free. Both are
reported for the real data.

## 4. Numerical verification
`results/tables/e10_occupancy_vs_rate.csv`, four ground truths:
- **OCC** clones differ in occupancy, shared rate → V_π > 0, V_logk = 0
- **RATE** clones share occupancy, differ in rate → V_π = 0, V_logk > 0
- **FROZE** clones differ in occupancy with k ≈ 0 (stable-state inheritance) → V_π > 0, V_logk = 0, no decay
- **BOTH** → both positive

OCC and FROZE have *identical occupancy distributions* and differ only in timescale, so the Milestone 2 plateau
is by construction the same for both. Separating them is the direct test of U1 and the reason the per-clone
estimator exists. A non-stationary arm applies a common logit shift to test U3's rescaling prediction.

**Results** (`results/tables/e10_focus.csv`, 20 clones, 10 replicates each):

| Mechanism | true V_π | est V_π | true V_log k | est V_log k | LRT shared π rejected | LRT shared k rejected |
|---|---|---|---|---|---|---|
| OCC occupancy varies | 0.033 | **0.030** | 0 | **0.000** | 100% | **0%** |
| RATE rate varies | 0 | 0.005 | 1.014 | 0.109 | 30% | **100%** |
| FROZE occupancy varies, k ≈ 0 | 0.033 | **0.178** | 0 | **3.14** | 100% | **0%** |
| BOTH | 0.027 | 0.029 | 0.947 | 0.172 | 100% | 90% |

**What works.** OCC and RATE are cleanly separated, which is the result the estimator exists for. The
likelihood-ratio test for a shared rate is the reliable discriminator: **0% false positives** under both
mechanisms whose true V_log k is zero (OCC and FROZE), 100% power under RATE and 90% under BOTH.

**What fails, and it is the predicted failure.** Under FROZE — clones frozen in different states, k ≈ 0.01, far
below the identifiability window — the noise-corrected point estimates break badly: V̂_π is five-fold too large
(0.178 against 0.033) and V̂_log k is spuriously enormous (3.14 against 0). The LRT nevertheless reports no
evidence for clone-specific rates, correctly. So:

> **When the rate falls outside the window of U4, the variance-component point estimates are not trustworthy,
> but the likelihood-ratio tests remain correctly calibrated.** Conclusions should therefore lead with the tests
> and treat the point estimates as descriptive. This is what the real-data analysis does.

**One caveat against the estimator.** Under RATE, V̂_π is slightly positive when the truth is zero (0.005) and the
LRT for shared π rejects in 30% of replicates rather than 5%. Rate heterogeneity therefore leaks a little into
the occupancy test. The real-data occupancy signal (p < 1e-10, V̂_π ≈ 0.03–0.05) is far larger than this leak,
but it is a known anti-conservatism of the method.

**Both components are under-estimated when true.** V_log k recovers 11–18% of the truth under RATE and BOTH at
this clone count, so the method is conservative about magnitudes even when it correctly detects the mechanism.

## 5. Cross-technology limits (the binding constraint)
From the published trees of two recording systems, the identifiable window for k spans only **0.09–0.24 decades**
in every dataset examined — macsGESTALT and Cassiopeia alike (see `milestone3_sampling_design.md`, table 2).
Switching rates are therefore not identifiable from either technology unless the true rate happens to fall inside
a window a factor of about 1.5 wide. This is a property of the recorders, not of the estimator, and it caps what
any method — ours, PATH, PASTRI or SCOUT — can extract about dynamics from these data.

Consequently, for the real data the defensible target is **V_π, clone-associated occupancy**, and claims about
V_logk must be presented as bounded by the window rather than as estimates.
