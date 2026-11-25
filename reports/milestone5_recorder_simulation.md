# Milestone 5: the recorder-realistic generative benchmark

Code: `src/branching/recorder.py`. Experiment: `experiments/reconstruction/run_e17.py`.

## 1. The observation hierarchy, implemented

| Level | What it is | How it is produced |
|---|---|---|
| **T0** | True genealogy: every cell, including those never sampled | `simulate_genealogy` — synchronous divisions over `n_gen` generations with per-cell death, optional clone-specific division probability, and an optional founder bottleneck that keeps only a few lineages at a chosen generation |
| **T1** | Sampled genealogy: ancestry among sampled cells, known without error | `induced_subtree` — the induced subtree on sampled leaves, unary chains collapsed, **retaining the true branch length in divisions** on every edge |
| **T2** | Recorded information: what the recorder actually writes | `record_barcodes` — cumulative per-site editing, homoplasic outcome distribution, dropout, readout missingness |
| **T3** | Reconstructed genealogy | `reconstruct` — `greedy`, `nj` or `upgma` on the T2 matrix |

Keeping true branch lengths at T1 is what makes the decomposition possible: comparing a fit on T1-with-true-lengths
against the same topology with unit edges isolates **branch-length miscalibration**, and comparing that against
T3 isolates **topology error**.

## 2. Biological processes represented
Cell division and death; clone-specific growth (via a per-founder division probability); metastatic founder
bottlenecks; unequal lineage expansion (an emergent consequence of death and bottlenecks); incomplete cell
sampling; and two-state phenotypic switching evolving along the true genealogy at rate `k` per generation.

Both **balanced** (death only) and **highly imbalanced** (bottleneck of 4 lineages at generation 4) genealogies
are generated, because imbalance is the property that drives the star-likeness we are testing.

## 3. Recorder processes represented
Cumulative editing of initially unedited sites; **saturation** (an edited site cannot be re-edited, so
information stops accruing); **homoplasy** (outcomes are drawn from a skewed, Zipf-like distribution so the same
indel recurs independently in unrelated lineages); **dropout** (a site silenced in a lineage); **readout
missingness**; and finite **barcode capacity** (number of sites).

## 4. Parameter provenance, and what is *not* claimed

Two settings are parameterised to resemble the two published systems. They are **plausible stand-ins, not fitted
reconstructions of either chemistry**, and no claim is made that they reproduce the published recorders
quantitatively.

| Parameter | `gestalt_like` | `cassiopeia_like` | Basis |
|---|---|---|---|
| recording sites | 9 | 25 | macsGESTALT uses a ~10-target evolving barcode; the Cassiopeia/A549 system carries many more integrated target sites. Orders of magnitude taken from the published designs; exact per-cell counts vary and are not fixed here. |
| edit probability per site per generation | 0.30 | 0.12 | Chosen so that each system's realised saturation over 10 generations is in a plausible range (see below). **Not measured from the data**; this is the least defensible parameter and is varied in sensitivity analysis. |
| distinct indel outcomes | 12 | 40 | Indel spectra are heavily skewed in both systems; a larger outcome space lowers homoplasy. Relative ordering is defensible, absolute values are not. |
| dropout | 0.01 | 0.005 | Nominal. |
| readout missingness | 0.12 | 0.20 | Both systems report substantial missing-allele rates; 10–25% is the range commonly quoted. |

**Realised recorder state in the simulation** (median over replicates), which is the quantity that actually
matters for inference:

| setting | saturation | missing | informative characters |
|---|---|---|---|
| `gestalt_like` | **0.97** | 0.19 | 42 |
| `cassiopeia_like` | 0.73 | 0.24 | 190 |

The `gestalt_like` setting saturates almost completely within ten generations — nearly every site is edited, so
late divisions leave no record. This is a direct, quantitative statement of why a low-capacity recorder cannot
resolve recent coalescences, and it is the mechanism behind the capacity law's ceiling.

## 5. Reconstruction algorithms
The published pipelines are not installable in this environment, so three standard algorithms are implemented
and used as stand-ins:

- **`greedy`** — recursive splitting on the most frequent non-ancestral character state. This is the algorithm
  Cassiopeia's greedy solver uses, reimplemented from its description; it is the closest available analogue to
  the published method.
- **`nj`** — neighbour joining on normalised Hamming distance over jointly observed sites.
- **`upgma`** — average linkage on the same distance.

Using three algorithms rather than one is deliberate: it supplies the ensemble for uncertainty-aware inference
and shows which conclusions are algorithm-specific.

## 6. Ground truth retained
Every replicate stores the complete genealogy, the sampled genealogy with true branch lengths, the true state of
every cell, the recorded allele matrix, and each reconstructed tree — so any later question about where
information was lost can be answered without re-simulating.

## 7. Honest limits
- Divisions are synchronous; real tumours are not.
- Editing is independent across sites, so inter-site dependence and large deletions spanning sites are absent.
  Cassiopeia-type systems do suffer multi-site dropout, which this under-models.
- The edit rate is the weakest parameter and is not derived from data.
- Reconstruction algorithms are reimplementations, not the authors' code.

Consequently this benchmark supports **relative** statements — how much information is lost between levels, and
how that depends on capacity and imbalance — and not absolute claims about either published system's recorder.
