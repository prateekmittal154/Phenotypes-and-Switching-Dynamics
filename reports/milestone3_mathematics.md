# Milestone 3: mathematics

Full statements and proofs: `theory/occupancy_vs_switching.md`. Estimator:
`src/lineage_inference/perclone.py`. Numerical verification: E10, E11.

## 1. The audit: three Milestone 2 claims were too strong

The brief asked whether the plateau identifies clone-associated occupancy rather than switching propensity. It
does. The algebra of the Milestone 2 decomposition is correct; the interpretation was not.

| Milestone 2 claim | Verdict | Correction |
|---|---|---|
| The plateau B = Var_c(π_c) measures **heritable switching propensity** | **Wrong** | B is invariant to every k_c. It measures **clone-associated equilibrium occupancy** (U1). Two populations with k differing 1000-fold have the identical plateau 0.0900. |
| B is **immune to stable-state selection** | **Conditionally wrong** | True only if clones share occupancy. Clones frozen in *different* states (k → 0, π_c ∈ {0,1}) give the maximal plateau, 0.25. The Milestone 2 simulation gave all clones the same rates, so it never tested this (U2). |
| B is **immune to environmental induction** | **Weaker than stated** | A common induction shift cannot *create* a plateau from zero — that part survives — but it does *rescale* an existing one: π_c ∈ {0.2, 0.8} gives B = 0.0900, and after a logit shift of 1 it gives 0.0653 (U3). Cross-compartment comparisons must use the logit scale. |

What this does **not** overturn: the Milestone 2 plateau estimates remain valid estimates of occupancy variance,
and the Simeonov finding that metastatic EMT enrichment is clonal composition is untouched — that analysis never
used the plateau.

## 2. The general model
Clone-specific forward and reverse rates (a_c, b_c), hence clone-specific occupancy π_c and timescale k_c;
non-stationary ancestral states; organ-specific shifts; unequal branch lengths; founder bottlenecks; missing and
non-random sampling. Reparameterised as (π_c, k_c) because the identifiability structure is clean in those
coordinates and badly coupled in (a_c, b_c).

## 3. The two results that matter

**U4 — joint identifiability, with an explicit window.** For a clone observed at stationarity with leaf pairs at
two or more distinct lineage distances, (π_c, k_c) is identified: the leaf marginal gives π_c, and the pair
covariance at distance L gives ρ = e^{−2k_cL} and hence k_c. Two distances over-identify it and make the model
testable. The rate is identifiable only inside

  **1 / (2 L_max) ≲ k_c ≲ 1 / (2 L_min),**

outside which the curve is flat (k indistinguishable from 0) or saturated (indistinguishable from ∞). This is the
per-clone sharpening of Milestone 2's Proposition T3, and it explains why that milestone's decay rate hit its
bound in 7 of 8 compartments: the real trees do not span the required window.

With per-clone estimates the two mechanisms finally separate:
- **V_π = Var_c(π_c)** — clone-associated occupancy, i.e. the Milestone 2 plateau under its correct name;
- **V_logk = Var_c(log k_c)** — clone-associated switching dynamics, which is what "heritable switching
  propensity" should mean and which Milestone 2 never estimated.

**U5 — the optimal lineage distance.** The Fisher information one within-clone binary pair at distance L carries
about log k is I(L) = I_ρ(ρ)(2kLρ)², ρ = e^{−2kL}. It vanishes at both ends and peaks at an interior optimum:

  **k L\* ≈ 0.40** (verified at k = 0.1, 0.5, 2.0 and π = 0.3, 0.5; the simplified L²e^{−4kL} form gives 0.5).

Information about π_c instead comes from the leaf marginal and scales with the number of **cells**. Hence the
design rule: **occupancy is bought with cells, switching rate is bought with pairs separated by about one
switching time.** Cells sequenced at an uninformative depth buy almost nothing about dynamics.

## 4. Honest placement
U1–U3 are corrections to my own earlier over-reach, not discoveries. U4 is elementary once the estimand is posed
in (π, k) coordinates; its contribution is the explicit identifiability window and the warning it issues about
real data. U5 is a standard Fisher-information calculation; its value is that it yields a concrete, falsifiable
design rule where the field's instinct is to sequence more cells. None of this is a new theorem, and none is
presented as one. PATH targets a population-average rate; PASTRI targets stage-specific rates; SCOUT is per gene.
None separates V_π from V_logk across clones, which is the gap these results fill.
