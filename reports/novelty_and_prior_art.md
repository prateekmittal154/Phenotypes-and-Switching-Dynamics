# Novelty and prior art (Milestone 1, Task 1)

**Search.**
- Consensus semantic search, 11 queries: selection vs plasticity, snapshot identifiability, switching-rate inference,
  sensitivity models, barcoding, trajectory-inference theory, unbalanced OT, Price equation, fitness and switching.
- PubMed metadata for the human atlases.
- GEO for data.
- This is a scoping audit, not a systematic review.
- Items marked † are cited from domain knowledge and were **not re-retrieved** in this run.

## 1. Directly competing work

| Work | Exact question | Measurements | Key assumptions | Identifiable / estimated | Non-identifiability acknowledged? | Bounds? Sharp? | Assumptions testable? | Validation | What remains open |
|---|---|---|---|---|---|---|---|---|---|
| Weinreb et al., *PNAS* 2018 | Which dynamics can be inferred from one snapshot | scRNA snapshot | Population balance, steady state, known growth | Dynamics given growth and density; growth vs transitions not separable | **Yes, central** | No bounds | Partly | Haematopoiesis fate assays | No primary/metastasis setting; no quantitative bound on the transition share |
| Schiebinger et al. 2019 (WOT)† and the unbalanced/dynamic OT family: TIGON (Sha 2023), Chizat 2022, SUDO (Ying 2026), stVCR (Peng 2026) | Couplings between timepoints with growth | Time-series snapshots | Entropic or WFR penalties; growth prior or penalty | Point estimate chosen by the penalty | Implicitly. The growth prior or penalty is the identification device. Ying 2026 shows concave penalties *separate* growth from transport. | No | No | Simulations, development | The penalty decides selection vs transport. Nobody reports the **range** over penalties. |
| Lavenant et al., *AAP* 2024; Guan et al. 2025 | Consistency and identifiability of SDE inference from marginals | Many timepoints | Gradient drift, known (Lavenant) or jointly inferred (Guan) diffusion; **no growth** | Drift (+ diffusion) iff out of equilibrium | Yes | No | Partly | Synthetic | Growth and selection excluded by assumption |
| Gunnarsson, Leder & Zhang, *J Theor Biol* 2023 | Switching and proliferation rates | Sorted-population fractions or counts over time | Multitype branching | Switching rates from fractions; net growth needs counts | Yes, practical | Likelihood CIs | Yes, in vitro | Cell lines | Needs sorted time courses; not available for metastases |
| Browning et al., *PLoS CB* 2025 | Discrete vs continuous adaptation | Low-cell-count proliferation | Stochastic phenotype PDE | Heterogeneity practically non-identifiable | Yes | No | Yes | In vitro | In vitro only |
| Piras et al., *Cell Rep Methods* 2026 | Active transitions from abundances alone | Cell-type counts over time | Multitype branching | Transition probabilities | Not emphasised | No | Partly | In vitro | Needs time series |
| Whiting et al., *Nat Commun* 2025 | Pre-existing vs switching resistance | Barcode lineage + population size over passages | Resistance-phenotype model | Pre-existence vs switching per cell line | Partly | Posterior | Yes (scDNA/scRNA follow-up) | 5-FU colorectal lines | Needs serial barcode passages; in vitro |
| Schiffman et al. (PATH), *Nat Genet* 2024 | Heritability vs plasticity of states | Phylogeny + state per cell | Markov states on a tree | Heritability, transition rates (point) | Partly | No | Partly | KPC mouse, GBM, B-ALL | No primary-vs-metastasis selection bound |
| Simeonov et al. (macsGESTALT), *Cancer Cell* 2021 | Metastatic potential by EMT state | Lineage + scRNA, primary + metastases (mouse) | Phylogeny ⇒ selection | "Late-hybrid EMT states aggressively selected" | No formal treatment | No | — | — | Is the selection claim **forced**, or are within-clone shifts compatible with induction? Testable with P7. |
| Quinn 2021; Yang 2022; Jones 2026 (KP-Tracer); Park 2026 (DNA Typewriter, bioRxiv) | Rates, routes, plasticity and pre-existing metastatic state | Evolving lineage + scRNA | Tree reconstruction | Dissemination topology; Park: founder state *before* transplant | Implicit | No | — | Functional | Park 2026 is the first design that records the founder state; it effectively observes γ(x, y) |
| Kim et al., *Cell* 2018 | Chemoresistance: selection of genotypes vs reprogramming | scDNA + scRNA, longitudinal | Genotype = clone | "Genotypes selected, transcriptomes reprogrammed" | No | No | — | — | Same question in treatment, without bounds |
| Goyal et al. (FateMap), *Nature* 2023; Iyer et al., *PLoS CB* 2025 | Pre-existing vs drug-induced fate | Barcodes + scRNA; time-lapse | Sister-cell concordance | Pre-determined fraction | Partly | No | Yes | In vitro | In vitro |
| Tagore et al., *Nat Med* 2025 (atlas audited here) | Brain-metastasis biology | snRNA + TCR + pooled lp-WGS + Slide-seq; **unpaired** primaries vs brain metastases | — | CIN-high state "pre-exists in primaries but enriched in brain metastases" | No | No | — | mIF in a paired cohort | By P4 this is compatible with selection but does not exclude induction. It is also unpaired, so the per-patient estimand is undefined. |
| Tan 2006†; Dorn, Guo & Kallus, *JASA* 2024; Masten et al. 2025; Zhang et al., *Biometrika* 2026 | Causal effects under bounded unobserved selection (marginal sensitivity model) | Observational data | Odds ratio of selection ≤ Λ | **Sharp** bounds on ATE/QTE/DTE | Yes | **Yes, sharp** | No (Λ is a sensitivity parameter) | Theory + data | **The same mathematical structure.** A bounded density ratio of selection gives sharp bounds; Λ is unverifiable. |
| Huber & Strassen 1973†; DeRobertis & Hartigan 1981†; hockey-stick divergence (Polyanskiy, Poor & Verdú 2010†) | Minimax tests; density-ratio classes | — | Density band {l ≤ π ≤ u} | TV geometry of bands | — | — | — | — | The **box-projection step of P2 is classical** in this family |
| Lindvall 1992† (maximal coupling) | min P(X ≠ Y) over couplings = TV | — | — | — | — | Sharp | — | — | Step 2 of P2 is classical |
| Luria & Delbrück 1943†; Shaffer 2017† | Pre-existing vs induced mutations/states | Variance across parallel cultures | Independent cultures | Mechanism via jackpot variance | — | — | Yes | Classic | P10: it does **not** transfer to metastases from one shared primary |
| Price 1970† | Selection vs transmission decomposition | — | Mean-trait change | Exact covariance decomposition | — | — | — | — | Mean-level only. τ is a distributional transmission analogue. |

## 2. Candidate contributions (≥ 3, compared before selection)

### Candidate A: sharp lower bound on induced transitions under bounded selection (P2–P4)
1. **Closest existing work.**
   - The marginal-sensitivity-model sharp bounds: Tan 2006†; Dorn et al. 2024.
   - The maximal-coupling theorem.
   - The density-band TV geometry: Huber–Strassen†.
2. **Why these do not solve the problem.** None is posed for founder-to-metastasis state transitions, and none gives
   the τ-curve.
3. **Improvement.** A closed-form, sharp τ_min(Λ) in hockey-stick form, interpolating TV (Λ = 1) and support mismatch
   (Λ = ∞).
4. **Assumptions.**
   - A bound Λ on founder-level selection.
   - A mean-field metastasis.
   - P measured at dissemination.
   - Representative sampling.
5. **Counterexamples that weaken it.** Yes:
   - founder bottleneck (P8);
   - regional sampling (P9);
   - plug-in undersampling (P5a).
   These weaken it, but do not invalidate the theorem.
6. **Data to test it.** No human data (G6). Mouse lineage data could test it (G7).
7. **What a meaningful negative result would look like.** The bound is ≈ 0 for every plausible Λ in realistic settings.
   E2 shows this for modest-selection families.

**Verdict:** mathematically correct, but a **direct composition of two classical results**. The novelty would be in the
application, not the mathematics.

### Candidate B: clone-constrained partial identification (P7)
1. **Closest existing work.** PATH (Schiffman 2024), Kim 2018, FateMap, Whiting 2025.
2. **Why these do not solve the problem.** They give point estimates or qualitative splits under a model, not a sharp
   set under within-clone selection bounds.
3. **Improvement.** Exact decomposition τ_min = Σ_c Q(c) τ_min,c. Without clone labels, the same assumptions are
   vacuous (= SM). With clone labels there is a Simpson reversal (C2), and state-pure clones add nothing (C3).
4. **Assumptions.** Exact heritable labels independent of expression; within-clone Λ_in; transitions preserve the
   clone.
5. **Counterexample.** Clone-label error biases the bound *upward* in principle, and the E4-S5 simulation must be read
   with care. inferCNV-type "clones" derived from the same RNA break independence.
6. **Data to test it.** **Yes, in mice:** GSE173958 (macsGESTALT) and GSE161363 (Quinn). **No, in humans.**
7. **What a meaningful negative result would look like.** In real lineage data, within-clone primary and metastatic
   state distributions are reconciled by Λ_in ≤ 5. The clone data then cannot exclude pure selection, and published
   "selection" claims are compatible with induction at the same strength.

**Verdict:** the most *useful* application. Mathematically elementary given A.

### Candidate C: minimal additional measurements (P4, P8–P10, P12)
1. **Closest existing work.** Weinreb 2018 (experimental-design remarks), the Luria–Delbrück fluctuation test, Park
   2026 (founder recording).
2. **Why these do not solve the problem.** No quantitative comparison across designs for the metastasis question. The
   failure of the fluctuation logic for a shared primary (P10) and the multiplication of Λ by regional bias (P9) are
   not stated in the retrieved literature.
3. **Improvement.**
   - The design law n_P ≈ ln(1/α′)(Λ(1 − q) + q)/q.
   - "More cells from one region make it worse."
   - Founder-level induction ≡ selection across any number of metastases.
   - Only founder-linked recording identifies τ without Λ.
4. **Assumptions.** The same as A, plus the founder-sampling model.
5. **Counterexample.** A memory effect of the pre-switch state on expansion breaks P10. That is a hidden-state
   enlargement, which the P10 proof covers.
6. **Data to test it.** Partly (mouse lineage data). The design law is checked by simulation (E5).
7. **What a meaningful negative result would look like.** None of the cheap measurements helps. That *is* the expected
   result for several designs, and it is publishable as a design principle.

**Verdict:** the most decision-relevant. Each statement is elementary. The **package** is not found in the retrieved
literature.

### Candidate D: continuous-time selection–plasticity identifiability (P11)
1. **Closest existing work.** Gunnarsson et al. 2023; Guan et al. 2025; Fassoni 2026; Qian 2014.
2. **Why these do not solve the problem.** They essentially do solve it: ≥ 3 out-of-equilibrium timepoints identify
   (a, c, Δb), and equilibrium sampling loses the split.
3. **Improvement.** A finite-state restatement only.
4. **Assumptions.** Homogeneous rates and noise-free fractions.
5. **Counterexample.** C4 and C5 confirm the known limits.
6. **Data to test it.** No longitudinal human metastasis single-cell data were found.
7. **What a meaningful negative result would look like.** Already known.

**Verdict: not new. Dropped.**

## 3. Selection
- **Primary direction: A + C, with B as the empirical vehicle.** Report the Λ-sensitivity curve and its design
  consequences (P2–P5, P7–P10, P12). The test bed is clone-resolved mouse lineage data.
- **Novelty claim allowed.** At most "a sharp, computable sensitivity analysis and a design calculus for the
  selection-vs-induction question in metastasis, built from classical tools".
- **Novelty claim not allowed.** "A new theorem". P2 is the maximal coupling theorem composed with a density-band TV
  projection. Its partial-identification logic is that of the marginal sensitivity model.

Gate G1 (frozen) is therefore **not passed** as a mathematical-originality gate. See `milestone1_decision.md`.

## References retrieved in this run (Consensus/PubMed)
- [Weinreb 2018](https://consensus.app/papers/details/w:ad05f9a4-2da4-4e23-8c99-2720bf006288/)
- [Sha 2023, TIGON](https://consensus.app/papers/details/w:aeab8353-197f-476b-b9b8-baf022d37c4f/)
- [Zhou 2021, MuTrans](https://consensus.app/papers/details/w:e225d194-d451-48fd-9d57-371d06da7247/)
- [Gunnarsson 2023](https://consensus.app/papers/details/w:3e7b74f2-ef6d-4a3f-acb0-2032b4598ab8/)
- [Browning 2025](https://consensus.app/papers/details/w:dd5b475d-6c1b-467e-b8ed-62865c03656b/)
- [Piras 2026](https://consensus.app/papers/details/w:33bd21d4-812b-4144-b4a4-73074a6ceee8/)
- [Simeonov 2021](https://consensus.app/papers/details/w:367976ef-20f8-46e7-9f4a-9dab50ec46a7/)
- [Schiffman 2024, PATH](https://consensus.app/papers/details/w:8d9e5840-d40b-4aa9-9a6b-68538e0a558c/)
- [Park 2026](https://consensus.app/papers/details/w:f84d1046-2309-455d-aaf8-dc23ac663df1/)
- [Masten 2025](https://consensus.app/papers/details/w:fa591b8a-8a60-42e6-828b-647f0af6ac3f/)
- [Zhang 2026](https://consensus.app/papers/details/w:0d0c3a75-7129-44ed-9291-d087771f2ac0/)
- [Dorn 2024](https://consensus.app/papers/details/w:e22aad9a-7be0-4e4d-8370-9a810f9bee54/)
- [Goyal 2023](https://consensus.app/papers/details/w:a5e4f49e-6faf-4e39-b06f-1543cc3b0f5b/)
- [Kim 2018](https://consensus.app/papers/details/w:2df20a89-f58c-472a-922c-c9be21d8a1a6/)
- [Lavenant 2024](https://consensus.app/papers/details/w:9a5793f7-7963-4a64-a5b8-800a780c3244/)
- [Chizat 2022](https://consensus.app/papers/details/w:e6027847-82c2-43e3-817f-51d0cb7d6ff5/)
- [Guan 2025](https://consensus.app/papers/details/w:486666de-12b4-45dd-9637-e13c0fa822d0/)
- [Fassoni 2026](https://consensus.app/papers/details/w:088b4f33-5d64-43cd-913d-906863de5cf4/)
- [Qian 2014](https://consensus.app/papers/details/w:5b28f115-c3ca-45b4-8547-ea28593b0300/)
- [Ying 2026, SUDO](https://consensus.app/papers/details/w:f7a43d92-9969-4cb8-bb9b-c06ebc0e4fed/)
- [Peng 2026, stVCR](https://consensus.app/papers/details/w:c1ed9a02-eb88-4e83-a098-732ffbbc3a6d/)
- [Stuart 2025, SCOUT](https://consensus.app/papers/details/w:eb9e359e-fd8a-4859-8d16-c6fda0723f1f/)
- [Iyer 2025](https://consensus.app/papers/details/w:2d8dc1b0-64f4-48a6-b307-6e39d6af8ba7/)
- [Whiting 2025](https://consensus.app/papers/details/w:046c999f-9a77-4982-82d9-aea1e0707b44/)
- [Yang 2022](https://consensus.app/papers/details/w:f3ddf678-8062-4322-b5b4-83115fe1f838/)
- [Quinn 2021](https://consensus.app/papers/details/w:f8183f00-d474-4edf-9628-5d28237cfd68/)
- [Jones 2026](https://consensus.app/papers/details/w:43efc914-a1b5-4fea-8f5f-fdce3dc4ca02/)
- Tagore et al., *Nat Med* 2025, [doi:10.1038/s41591-025-03530-z](https://doi.org/10.1038/s41591-025-03530-z) (PubMed)
- Burns et al., *JCI* 2026, [doi:10.1172/JCI194708](https://doi.org/10.1172/JCI194708) (PubMed)
