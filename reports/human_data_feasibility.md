# Human-data feasibility (Milestone 1, Task 7)

**Sources.** GEO SOFT and family files were downloaded on 2026-11-23 into `data/provenance/`, and parsed. No expression
matrix was downloaded: the audit needed metadata only. The matrix is `data/manifests/dataset_suitability_matrix.csv`.

## NSCLC brain-metastasis atlas (GSE223503 SuperSeries; Tagore et al., *Nat Med* 2025)
**Components (verified):**
- GSE223499: snRNA-seq, 42 samples.
- GSE223500: TCR, 42 samples.
- GSE223501: Slide-seqV2, 14 samples (8 primary tumours, 6 brain metastases).
- GSE223502: low-pass WGS, 38 samples.

**Pairing: none.**
- Primary-tumour patient IDs, e.g. N254, N586, KRAS_10–13, STK_15/20/21/22: 11 in snRNA.
- Brain-metastasis IDs, e.g. PA001–PA141, KRAS_6/7/8/17, STK_1/3/5/14/18: 31 in snRNA.
- **No ID overlap** in any component.
- The paper's paired validation used multiplexed immunofluorescence in an *independent* cohort. It is not a public
  single-cell resource.

**Genotype.**
- The lp-WGS is per sample (pooled), not per cell.
- Per-cell "clones" would have to be inferred from the snRNA itself (inferCNV). They are **not independent** of the
  phenotype measurement, and P7 requires independence.

**Treatment.** Treatment-naive by design (series summary). This helps, but it is moot without pairing.

**Access.** Processed counts are public. Raw data are declared for dbGaP.

**Metadata quality.** Four GSE223499 samples named KRASmut_BM_6–9 carry "metastatic site: PT". This is a naming
inconsistency to resolve if the atlas is ever used.

**What it can support.**
- Population-level, unpaired comparisons only. τ is defined per patient, through the per-patient founder→metastasis
  coupling, so pooling unpaired primaries and metastases changes the estimand to "metastasis-cohort vs primary-cohort
  state frequencies". Selection *of patients*, i.e. which patients get resected brain metastases, then enters as
  another unbounded selection term.
- The paper's claim that "the CIN-high subpopulation pre-exists in primaries but is strongly enriched in brain
  metastases" is, by P4, **compatible with** pure selection and **does not exclude** induction. It cannot be checked
  at patient level with these data.

## GSE248830
- Breast and lung cancer with paired brain metastases.
- FFPE **bulk microarray**: 22 LUAD + 22 BC samples.
- No cell-state distributions.

## GSE271259 (Burns et al., *JCI* 2026)
- Bulk RNA-seq: 5 author-declared LUAD/brain-metastasis pairs, plus MET-amplified vs non-amplified brain metastases.
- SRA not public (consistent with the earlier metastasis-pathways audit).
- Not single-cell.

## Gate G6: **FAIL**
- No audited human dataset gives patient-matched primary and metastatic malignant single-cell state distributions for
  ≥ 5 patients.
- No source calibrates Λ for transcriptional states in human metastasis.
- Human application is deferred. This was anticipated: the SAP prioritises synthetic and experimental ground truth.

## Gate G7: **PASS** (existence and access verified; content not yet inspected)
- **GSE173958, Simeonov et al. 2021 (macsGESTALT, KPC PDAC).**
  - 2 mice, each with primary tumour, surgical site, peritoneal metastases, liver, lung and CTCs.
  - scRNA, plus static and evolving CRISPR barcodes in per-sample `*.stats.txt.gz` lineage files.
  - Public (1.67 GB).
  - The direct test bed for P7: the clone labels come from DNA barcodes, independent of expression.
- **GSE161363, Quinn et al. 2021 (A549 xenograft).**
  - Allele tables, character matrices, trees and scRNA for 4 lineage libraries.
  - Public (820 MB).

**Caveats for both.**
- The primary is sampled at the endpoint, not at dissemination: the P measurement-time mismatch [A].
- Founder states are not recorded.
- Park et al. 2026 (DNA Typewriter) would record founder states; its data availability was not verified.
