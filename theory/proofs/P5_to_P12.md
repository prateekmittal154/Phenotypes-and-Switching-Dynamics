# Proofs of P5–P12

## P5 (finite samples)
**(a) The plug-in bound is invalid.**
- Suppose P̂(s) = 0. Then for every m, A(m) = Σ(Q̂ − ΛmP̂)₊ ≥ Q̂(s), so τ̂_min ≥ Q̂(s).
- P(P̂(s) = 0) = (1 − p)^{n_P}.
- Under pure selection τ = 0, so the plug-in reports τ ≥ Q̂(s) > 0 with at least this probability, for every Λ,
  including the true one. ∎

**(b) The robust bound is valid.**
- Let E be the event that P ∈ box_P and Q ∈ box_Q. Each box is a Bonferroni-simultaneous Clopper–Pearson box at level
  α/2, so P(E) ≥ 1 − α.
- On E, τ_L = min over the boxes of τ_min(Λ; P′, Q′) ≤ τ_min(Λ; P, Q) ≤ τ, for every mechanism in M(Λ).
- **LP form.**
  - The constraint mP′ ≤ π ≤ ΛmP′ is bilinear in (m, P′). Substitute u = mP′.
  - P′ ∈ [L, U] with ΣP′ = 1 is equivalent to mL ≤ u ≤ mU with Σu = m.
  - The constraints become u ≤ π ≤ Λu and Σπ = 1.
  - TV(π, Q′) = min Σs subject to s ≥ Q′ − π and s ≥ 0.
- Every constraint is linear in (π, u, m, Q′, s), so τ_L is the value of one LP with 4d + 1 variables. m > 0 holds
  automatically, because Σπ = 1. ∎

**Note on unseen states.** A metastatic state never seen in the primary has L = 0 and U = 1 − (α′/2)^{1/n_P} > 0. The
LP may therefore place primary mass up to U on it. This is how Mechanism U (rare unsampled states) is built into the
bound.

## P6 (measurement error and Lipschitz constants)
**1-Lipschitz in Q.**
- For Q, Q′: |A_Q(m) − A_{Q′}(m)| ≤ Σ(Q − Q′)₊ = TV(Q, Q′), and the same holds for B.
- max and min (over m) are 1-Lipschitz in the sup norm.
- Hence |τ_min(P, Q) − τ_min(P, Q′)| ≤ TV(Q, Q′).
- For Q′ = QM: TV(QM, Q) ≤ Σ_x Q(x)(1 − M(x, x)) ≤ ε.

**Lipschitz constant in P.**
- Fix r ∈ [m, Λm] and write a = ⟨P, r⟩, a′ = ⟨P′, r⟩, δ₁ = ‖P − P′‖₁.
- Because Σ(P − P′) = 0, |a − a′| = |⟨P − P′, r − (1 + Λ)m/2⟩| ≤ δ₁(Λ − 1)m/2.
- Then ‖Pr/a − P′r/a′‖₁ ≤ (‖(P − P′)r‖₁ + |a − a′|)/a ≤ δ₁(Λm + (Λ − 1)m/2)/m = δ₁(3Λ − 1)/2.
- Each π ∈ Π_Λ(P) therefore has a partner in Π_Λ(P′) within TV ≤ TV(P, P′)(3Λ − 1)/2. The infimum inherits this
  constant. ∎

## P7 (clone-resolved bound)
**Setting.** With w(c, x) = v(c) r_c(x), the selected joint law is π(c, x) ∝ P(c, x) v(c) r_c(x).
- Transitions preserve the clone, so Q(c) = π(c).
- Within clone c, π(· | c) ∈ Π_{Λin}(P(· | c)).
- Because v is unconstrained, any clone marginal with supp ⊆ supp P(c) is attainable:
  v(c) ∝ Q(c) / (P(c) E_{P(·|c)} r_c).

**Decomposition.**
- τ = Σ_c Q(c) τ_c, where τ_c is the within-clone off-diagonal mass.
- τ_c ≥ τ_min(P(· | c), Q(· | c); Λ_in) by P2.
- The per-clone minimisers are chosen independently, so the sum is attained. ∎

**(a) No clone labels.** The population analyst must allow every latent clone structure consistent with (P, Q). This
includes one state-pure clone per state, under which v realises any reweighting of P on supp P. Hence τ_min^pop = SM(P, Q).

**(b) Simpson reversal.** See `../counterexamples/C2.md`.

**(c) State-pure clones.** Each conditional law is a point mass, so τ_c = 0.

**(d) Clone-label error.** Label error is not covered by the bound; E4-S5 shows the effect numerically.

## P8 (founder bottleneck)
- Let x₁…x_{N_f} be iid ~ P (neutral seeding, Λ = 1) with K = I, and Q = Σ Z_i δ_{x_i} / ΣZ_i. Take a target Q\* with
  k ≤ N_f atoms in supp P.
- The event that the founders occupy exactly those k states, in some assignment, has positive probability.
- Given that assignment, the map (Z_i) ↦ composition is continuous and onto the relative interior of the face spanned
  by those states. The Z_i have a positive joint density.
- So every neighbourhood of Q\* (relative to the face) has positive probability, and the realised Q follows no fixed
  ratio constraint to P.
- **Consequence.** Take any decision rule that declares induction on a set of Q values with non-empty interior in such
  a face. One example is the set where one rare primary state dominates. That rule has false-attribution probability
  bounded away from 0, uniformly in n_Q and n_P.
- **Mean-field recovery.** If the clone-size law does not depend on state, exchangeability gives E[Q(x)] = π(x). ∎

## P9 (regional sampling)
- P_obs = Pv / ⟨P, v⟩, so π ∝ Pw = P_obs (w/v) · const.
- If w has ratio ≤ Λ and v has ratio ≤ ρ, then w/v has ratio ≤ Λρ.
- Conversely, any reweighting of P_obs with ratio ≤ Λρ is attainable by some (w, v) pair. So the analyst of P_obs
  faces M(Λρ).
- ρ is a property of the region sampled, not of the number of cells, so it does not shrink as n grows. ∎

## P10 (founder-level induction is invisible to fluctuation)
**Mechanism A: founder-level induction.**
- Each metastasis draws N_f (any law) founders with primary states x_i iid ~ π.
- Each founder switches to y_i ~ K(· | x_i), independently and heritably.
- Expansion depends on y_i only.
- The y_i are then iid ~ πK.

**Mechanism B: pure selection.** Founders drawn iid ~ πK, with the same expansion law.

**Equivalence.** The joint law of (y_i, Z_i) is the same in both mechanisms. Metastases are conditionally independent
given the shared primary in both, so every joint law of compositions across metastases coincides.

**Caveat.** The equivalence fails if expansion depends on the *pre-switch* state x_i, a memory effect. That is a hidden
heritable variable: enlarging the state space to (x, y) restores the equivalence at the enlarged level. ∎

## P11 (dynamics)
**Model.**
- Counts follow n′ = n(B + T), with B = diag(b) and T a switching generator.
- The fractions p = n/Σn obey p′ = pA − (p·A1)p, with A = B + T.
- For d = 2, x = p₂: x′ = a(1 − x) − cx + Δb x(1 − x), where Δb = b₂ − b₁.

**(a) Two snapshots.**
- Pure selection (a = c = 0) gives a logistic curve: logit x(t) = logit x₀ + Δb·t.
- Pure induction (Δb = c = 0) gives x(t) = 1 − (1 − x₀)e^{−at}.
- Each one-parameter family hits any x₁ ∈ (x₀, 1) at a chosen t.

**(b) The vector field identifies A up to A + cI.**
- Suppose A and A′ give the same field on the simplex. Let D = A − A′ and δ = D1.
- Then pD = (p·δ)p for all p ∈ Δ.
- With p = e_i: D_{i·} = δ_i e_i, so D is diagonal with diagonal δ.
- With p = (e_i + e_j)/2: pD = (δ_i e_i + δ_j e_j)/2 and (p·δ)p = ((δ_i + δ_j)/2)(e_i + e_j)/2. Comparing the i-th
  coordinates gives δ_i/2 = (δ_i + δ_j)/4, so δ_i = δ_j.
- So D = δI.
- For d = 2, x′ is a quadratic in x with coefficients (a, Δb − a − c, −Δb). This map is one-to-one, so a trajectory
  observed at ≥ 3 non-equilibrium times identifies (a, c, Δb) generically.

**(c) Stationarity.** At a stationary x\*, the single equation a(1 − x\*) − cx\* + Δb x\*(1 − x\*) = 0 has a
two-parameter solution set in (a, c, Δb). See C5. ∎

**Prior art.** Gunnarsson, Leder and Zhang (J Theor Biol 2023) show that fraction data identify switching rates but not
absolute net growth. Guan et al. (2025) prove that, for gradient-flow SDEs, drift and diffusion are identifiable iff the
process is observed out of equilibrium. P11 is the finite-state analogue and is **not claimed as new**.

## P12 (design law)
**Maximal selected frequency.**
- Under M(Λ), the largest π(s) given P(s) = p is attained at w = Λ on s and w = 1 elsewhere. It equals
  Λp/(Λp + 1 − p).
- With K = I (pure selection), Q(s) = π(s). So pure selection explains Q(s) = q iff p ≥ p\*(q, Λ) = q/(Λ(1 − q) + q).

**Required number of primary cells.**
- If s is unseen in n_P primary cells, the one-sided Clopper–Pearson upper bound at level α′ is 1 − α′^{1/n_P}.
- Requiring this upper bound to be < p\* gives n_P > ln(1/α′)/(−ln(1 − p\*)).
- This is a sufficient single-state condition. The LP τ_L uses every state and agrees with it in E5.

**Lineage recording.**
- Suppose each metastatic cell carries an error-free record of its founder's state.
- Then observing one cell with a recorded founder state ≠ its current state proves τ > 0.
- P(at least one such cell in n) = 1 − (1 − τ)^n.
- With recording error ε_r, the test becomes a binomial test against ε_r. ∎
