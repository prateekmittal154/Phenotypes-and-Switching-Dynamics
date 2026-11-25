# Milestone 3: what the plateau actually identifies (Propositions U1–U6)

This section **corrects three claims made in Milestone 2**. The algebra of Proposition T2 is unchanged and remains
correct; its *interpretation* was too strong. Proofs below; numerical verification in E10 and E11.

Throughout, clone c has forward rate a_c, reverse rate b_c, equilibrium occupancy π_c = a_c/(a_c + b_c) and
switching timescale k_c = a_c + b_c. The decomposition of T2 is C(L) = E_c[π_c(1 − π_c)e^{−2k_cL}] + Var_c(π_c).

## U1. The plateau measures occupancy, not switching propensity

**Statement.** B = Var_c(π_c) is invariant to the switching timescales {k_c}. It identifies clone-associated
*equilibrium state occupancy*. It carries no information about whether a clone reaches that occupancy by rapid
balanced switching or by never switching at all.

**Proof.** In the law-of-total-covariance step, B is the variance of the conditional mean E[X | c], which at
stationarity equals π_c for every c, whatever k_c is. Nothing else in B depends on the generator. ∎

**Counterexample.** Two populations with the same occupancies π_c ∈ {0.2, 0.8} but k = 0.01 in one and k = 10 in
the other have the **same plateau, 0.0900**, and the same A(0) = 0.1600, while their switching dynamics differ by
three orders of magnitude.

**Consequence.** The Milestone 2 phrase "heritable switching propensity" is wrong for B. The correct name is
**clone-associated state occupancy**. Genuine heritable switching *dynamics* is Var_c(k_c), which the Milestone 2
estimator never estimated, because it fitted a single shared k.

## U2. Immunity to stable-state selection holds only under a further assumption

**Statement.** The plateau is immune to within-clone state memory, which lives entirely in A(L). It is **not**
immune to between-clone state differences, regardless of their origin. In particular, if "stable-state
inheritance" means clones are frozen in *different* states (k_c → 0, π_c ∈ {0, 1}), then B = Var_c(π_c) is at its
maximum, not zero.

**Proof.** Immediate from U1: B depends only on the spread of π_c. ∎

**Why Milestone 2 reported immunity.** In that simulation the stable-state scenario gave every clone the *same*
rate parameters, hence the same π_c, so the true B was zero by construction and the estimator correctly returned
zero. That experiment tested robustness to founder sampling and drift, which is a real and useful result, but it
did not test a stable-state mechanism with heterogeneous frozen states. E10's FROZE arm supplies that test.

## U3. Immunity to environmental induction survives only in a weaker form

Let induction apply a clone-independent shift δ_o on the logit scale, π_{c,o} = expit(logit π_c + δ_o).

**(a) Induction cannot create a plateau.** If Var_c(π_c) = 0 then all π_c are equal, a common shift keeps them
equal, and Var_c(π_{c,o}) = 0. ∎

**(b) Induction does rescale an existing plateau.** Since dπ′/dπ = π′(1 − π′)/(π(1 − π)), the delta-method gives

  Var_c(π_{c,o}) ≈ Var_c(π_c) · [π̄′(1 − π̄′) / (π̄(1 − π̄))]².

Worked instance: π_c ∈ {0.2, 0.8} has B = 0.0900; after δ = 1 the occupancies are {0.405, 0.916} and
B = 0.0653; after δ = 2 they are {0.649, 0.967} and B = 0.0254.

**Consequence.** Milestone 2's unqualified statement "induction contributes nothing to B" is false. The defensible
statement is (a). Comparing B between compartments therefore requires working on the logit scale, or reporting
Var_c(logit π_c), which is invariant to a common shift.

## U4. Joint identifiability of occupancy and switching rate within a clone

**Statement.** For one clone observed at stationarity with leaf pairs at two or more distinct lineage distances,
(π_c, k_c) is identified, provided 0 < π_c < 1 and the observed distances satisfy the window condition below.
With one distance only, k_c is identified but not testable. The identifiable window is

  **1/(2L_max) ≲ k_c ≲ 1/(2L_min).**

**Proof.** The leaf marginal gives π_c directly. For a pair at distance L the covariance is π_c(1 − π_c)e^{−2k_cL},
and since π_c is known, ρ(L) = e^{−2k_cL} is observed, so k_c = −log ρ(L)/(2L). Two distances over-identify k_c and
make the two-state model testable. Outside the window the map is flat: if 2k_cL_max ≪ 1 then ρ ≈ 1 at every
observed distance and k_c is indistinguishable from 0; if 2k_cL_min ≫ 1 then ρ ≈ 0 everywhere and k_c is
indistinguishable from ∞. ∎

This is the per-clone refinement of Milestone 2's Proposition T3, and it is what the per-clone estimator in
`src/lineage_inference/perclone.py` implements, with profile-likelihood intervals and an explicit flag when k̂
reaches a bound.

**Variance components.** With per-clone estimates in hand, the two mechanisms separate:
- **V_π = Var_c(π_c)** — clone-associated occupancy (the Milestone 2 plateau, correctly named);
- **V_logk = Var_c(log k_c)** — clone-associated switching dynamics, the quantity that actually deserves the name
  "heritable switching propensity".
Both are corrected for estimation noise by subtracting the mean squared standard error.

## U5. The optimal lineage distance for learning a switching rate

**Statement.** The Fisher information that one within-clone binary pair at lineage distance L carries about log k_c
is

  I_{log k}(L) = I_ρ(ρ) · (2 k_c L ρ)²,  ρ = e^{−2k_cL},
  I_ρ(ρ) = (π q)² [1/p₁₁ + 1/p₀₀ + 2/p₀₁],  p₁₁ = π² + πqρ, p₀₀ = q² + πqρ, p₀₁ = πq(1 − ρ).

It vanishes as L → 0 and as L → ∞, and is maximised at an interior optimum. Dropping the weak ρ-dependence of
I_ρ, the maximiser of L²e^{−4kL} is **L\* = 1/(2k)**; numerically, with I_ρ retained, the optimum sits at
**k L\* ≈ 0.40**, verified at k = 0.1, 0.5 and 2.0.

By contrast the information about π_c comes from the leaf marginal and grows with the **number of cells**,
essentially independently of depth.

**Consequence, and the design rule.** *Occupancy is bought with cells; switching rate is bought with pairs
separated by about one switching time.* Sequencing more cells at an uninformative lineage depth buys almost
nothing about dynamics. This formalises the Milestone 2 observation that depth spread mattered more than cell
count, and it is the central design result of Milestone 3.

## U6. What this implies for the Milestone 2 estimator
The two-step estimator fitted a single shared k. Under clone-specific k_c it recovers a weighted average rate; by
U1 this does not bias the plateau, so the Milestone 2 plateau estimates stand as estimates of **occupancy
variance**. What they do not support is the interpretation placed on them. The per-clone estimator supersedes it
for any question about dynamics.

## Relation to existing methods
- **PATH** (Schiffman 2024) estimates a single phylogenetic autocorrelation and converts it to transition rates
  under a shared-rate model. It therefore targets a population-average k and does not separate V_π from V_logk.
- **PASTRI** (Yang 2026) uses pairs at varying phylogenetic distance to obtain *stage-specific* rates, i.e.
  heterogeneity in k across developmental time rather than across clones.
- **SCOUT** (Stuart 2025) fits Ornstein–Uhlenbeck processes per gene on lineage trees; the OU reversion rate is the
  continuous analogue of k, but the model is per gene and does not decompose across clones.
- U4's window condition and U5's optimal distance are, as far as the searches found, not stated for this estimand;
  both are elementary once posed, and are presented as such rather than as new theorems.
