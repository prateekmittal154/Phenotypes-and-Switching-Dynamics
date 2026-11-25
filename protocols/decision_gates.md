# Study D — Milestone 1 decision gates (frozen before any simulation output was inspected)

These rules are frozen before `experiments/` is run. Their SHA-256 checksum is in `checksums/FREEZE_SHA256_m1.txt`.
Any later change is logged in `protocols/deviations.md` and does not alter the decision.

The *central result* is the primary direction selected in `reports/novelty_and_prior_art.md` (§ Selection).

| Gate | Operational criterion | PASS | FAIL consequence |
|---|---|---|---|
| G1 Originality | No retrieved work states the central result, or a mathematically equivalent one, in a form directly usable for primary-vs-metastasis inference. | No equivalent found **and** the result is not a direct corollary of one named classical theorem | If it is a corollary of a classical result with only re-interpretation: **REFORMULATE** |
| G2 Mathematical validity | Complete written proof. The closed form agrees with an independent LP to ≤ 1e-7 on ≥ 1,000 random instances. | Both | **STOP** |
| G3 Nontriviality | Across the frozen random-instance family and Λ ∈ {2, 5, 10, 100}, τ_min(Λ) lies strictly between the support-mismatch bound and TV(P, Q). | Gap ≥ 0.05 on both sides in ≥ 50% of instances for at least two Λ values ≤ 10 | **REVISE** |
| G4 Computational feasibility | One evaluation at d = 1,000 states. | Closed form < 1 s; robust LP < 60 s | **REVISE** |
| G5 Robustness | Finite-sample lower confidence bound (α = 0.05) under S1, S3, S4, S5 and S6 with correctly specified Λ. | Empirical coverage ≥ 0.93 in every scenario and false induction attribution ≤ 0.05 under pure selection (S1, S4) | **REVISE** (finite-sample procedure) |
| G6 Biological relevance | At least one audited public human dataset gives **patient-matched** primary and metastatic malignant single-cell state distributions for ≥ 5 patients. Λ must also be calibratable from an external source. | Both | Human application deferred. Not by itself STOP. |
| G7 Independent validation | At least one identified public dataset (any species) gives lineage-linked primary/metastatic states, i.e. clone or founder labels measured independently of expression. | Exists and is accessible (verified) | **REVISE** |

## Decision rule
- **STOP** if G2 fails, or if G1 fails **and** G3 fails, or if G1 fails **and** both G6 and G7 fail.
- **REFORMULATE** if G1 fails but G2–G4 pass. The contribution is then an application or design principle, not a theorem.
- **GO** only if G1–G5 pass and at least one of G6 or G7 passes.
- **REVISE** otherwise.

A partial pass is not rounded up. Thresholds are not relaxed after results are seen.
