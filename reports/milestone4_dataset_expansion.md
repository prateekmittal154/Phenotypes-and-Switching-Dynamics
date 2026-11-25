# Milestone 4: dataset expansion and provenance

Inventory of every lineage-plus-phenotype resource considered, what was actually retrieved, and which analyses
each can support. Provenance files are in `data/provenance/`; downloaded file sizes and checksums are recorded
where files were retrieved.

## A. Quinn et al., GSE161363 (Cassiopeia, A549 lung adenocarcinoma xenograft) — **retrieved and analysed**

**Replication units established from metadata, not from clone counts.**

| Mouse | cells with lineage | lineage groups | anatomical sites | primary-site cells |
|---|---|---|---|---|
| 5k | 37,814 | 100 | 6 (LL, RE, RW, M1, M2, Liv) | 11,106 |
| 30k | 17,189 | 4 | 6 | 1,588 |
| 10k | 2,541 | 13 | 9 | 930 |
| 100k | 5,806 | 31 | 6 (LL, NE, NW, SE, SW, MS) | 4,430 |
| **total** | **63,350** | **148** | — | **18,054** |

**Four independent animals**, against two in the Simeonov system. Each mouse carries a barcode library of
different complexity (the 5k/10k/30k/100k labels), so the four are not technical replicates of one another.
Sites comprise lung lobes and quadrants (primary implantation), mediastinal nodes (M1, M2, MS), liver (Liv) and
blood.

**Paired expression and lineage.** Barcode naming differs per mouse: 5k, 30k and 100k match the matrix barcodes
directly, 10k requires the `<site>.` prefix to be stripped. Matching was done by whichever convention recovered
more cells, giving the counts above (89–97% of each mouse's metadata rows).

**Trees.** Three reconstructions are deposited per lineage group (`hybrid_priors.alleleThresh`, `nj`,
`nj_priors`). Only the authors' primary hybrid-with-priors tree is used; using all three triple-counts every
cell, which an early run did before the error was caught (deviation D20). Trees above 1,500 matched cells are
skipped and at most 30 trees per mouse are fitted, to bound the likelihood-ratio optimisation.

**Tree depth.** Cassiopeia trees are substantially deeper than macsGESTALT: median maximum coalescent depth
8–14 against 3–4. This is the main reason the Quinn system is worth analysing — it probes slower switching
rates and covers about half a decade more of rate space.

**Cell states derived within this system (deviation D19).** The Milestone 4 brief requires state definitions
that do not import the Simeonov EMT axis. Two were derived:
- *Unsupervised PC1* of the top 2,000 detected genes (mitochondrial and ribosomal genes excluded), fitted on
  primary-site cells only, with library size regressed out. **This axis was rejected as the primary definition**:
  it correlates with mitochondrial fraction at −0.39, +0.60, −0.14 and −0.67 across the four mice, so it is
  substantially a cell-quality axis. It is retained as a labelled sensitivity analysis.
- *Cell-cycle score* from a fixed external S/G2M gene set, binarised at the primary-site median. This is
  biologically interpretable, is not fitted to these data, and is not circular with lineage. It is the primary
  state definition, and it also serves as a useful positive control: proliferation is expected to be dynamic,
  so a competent method should find low state heritability and a fast switching rate.

Agreement between the two definitions is 36–43%, i.e. they are largely independent axes.

## B. Spatial KP lineage data (Jones et al. 2026) — **not retrieved**
No GEO entry was found by searches for spatial lineage tracing in Kras;Trp53 lung adenocarcinoma, or for
KP-Tracer with spatial coordinates. The published dataset may be deposited under an accession not indexed by
these terms, or in a non-GEO repository. **Missing resource:** a public accession carrying matched lineage
barcodes, expression and spatial coordinates. Without it the spatial question is not addressable, and no
spatial analysis is reported. No substitute was forced.

## C. Park et al., DNA Typewriter metastasis study — **not retrieved**
Targeted GEO searches (`DNA Typewriter`, `Typewriter lineage recording`, plus metastasis- and ENO1-specific
terms) and an SRA search returned no matching dataset. Three DNA Typewriter datasets do exist publicly
(GSE341627, a whole-mouse lineage reconstruction; GSE315827, the SciPhy framework; GSE291244, gastruloid
heterogeneity), **none of which is the lung-cancer metastasis study**. The preprint was located in Milestone 1
but its data availability was never established.

**Precisely what is missing:** the accession for the time-resolved DNA Typewriter recordings of orthotopically
transplanted lung cancer cells with matched pre-transplant transcriptomes. Until it is public, the claim that
timing information resolves parameters static trees cannot is **untested here**. What can be said from theory:
the identifiability window of Proposition U4 is set by the resolvable coalescent depths, and a recorder that
timestamps edits converts depth into elapsed time, which removes the clock assumption but does not by itself
widen the window — widening requires more resolvable levels (the capacity law in
`milestone4_design_optimization.md`).

## D. Other cohorts screened
GSE341627 (whole-mouse DNA Typewriter) and GSE315827 (SciPhy) carry lineage recordings but no matched cancer
cell-state phenotype of the kind this framework needs, and GSE291244 is a gastruloid fluctuation experiment.
None was analysed. Priority for future expansion is a recorder chemistry different from both GESTALT and
Cassiopeia, with more resolvable levels.

## Eligibility matrix

| Dataset | occupancy V_π | switching V_log k | composition decomposition | design/window audit | spatial |
|---|---|---|---|---|---|
| Simeonov GSE173958 (2 mice) | yes | **no** — rates above the covered range | yes | yes | no |
| Quinn GSE161363 (4 mice) | yes | tested, see `milestone4_biological_validation.md` | yes | yes | no |
| Jones 2026 spatial | — | — | — | — | not retrieved |
| Park DNA Typewriter | — | — | — | — | not retrieved |

## Honest summary
The expansion doubled the number of independent animals (2 → 6 across two systems) and added a second recorder
chemistry with deeper trees. It did **not** deliver the time-resolved or spatial data the brief hoped for,
because those are not publicly locatable, and that is recorded as a missing resource rather than worked around.
