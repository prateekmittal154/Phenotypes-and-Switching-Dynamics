# Research questions (Milestone 1)

**Q0, the central question.** Given primary-tumour and metastatic malignant-cell state measurements, which evolutionary
explanations (selection S, induction I, mixed M, unobserved heterogeneity U) are compatible with the data? Which can be
excluded, and under what explicit assumptions?

**Q1, the estimand.** τ is the fraction of metastatic cells whose state differs from the state of the primary-tumour
founder cell they descend from. It is a lineage-level quantity, not a count of transition events. τ = 0 is "pure
selection".

**Q2, exclusion of pure selection.** Under a bound Λ on relative selection (dissemination × survival × expansion) of
founder states, what is the sharp lower bound τ_min(Λ)? How does it compare with TV(P, Q) and with the
support-mismatch bound?

**Q3, exclusion of induction.** Can any cross-sectional (P, Q) data exclude induction?

**Q4, finite samples.** Is the plug-in τ_min(Λ) a valid lower bound when rare primary states are undersampled? If not,
what valid lower confidence bound replaces it?

**Q5, measurement design.** Which additional measurements make τ identifiable or reduce the bound's dependence on Λ?
Candidates: more cells, more regions, clone labels, a second timepoint, multiple metastases, founder-linked lineage
recording. Compare them under explicit assumptions.

**Q6, human data.** Does any audited public dataset support Q2–Q5 with patient-matched primary and metastatic
single-cell states?

## Out of scope for Milestone 1
- Continuous state spaces.
- Spatial kernels.
- Any classifier of primary vs metastatic identity.
