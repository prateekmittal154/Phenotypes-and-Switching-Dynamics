# Milestone 5: where switching information is lost

Experiment E17 (`experiments/reconstruction/run_e17.py`), 480 replicates: 2 recorders × 2 genealogy types ×
3 true rates × 40 replicates, 96 sampled cells per replicate. Theory: `theory/information_hierarchy.md`.

The same observed leaf states are fitted four ways, so the differences isolate one error source at a time:

| Level | What the analyst conditions on | Adds |
|---|---|---|
| **L1** | true sampled topology **with true branch lengths** | — (reference) |
| **L2** | true sampled topology, **unit edges** | branch-length miscalibration |
| **L3** | **reconstructed** topology, unit edges | topology error |
| **L3-ens** | likelihood averaged over three reconstructions | uncertainty integration |

## Result 1 — the dominant error is treating edit depth as a clock, not getting the topology wrong

Paired within-replicate log-ratio of $\hat k$ against the L1 fit (0 = no change):

| Level | median log-ratio | fold error | IQR | fraction overestimating |
|---|---|---|---|---|
| **L2** (branch lengths only) | **+0.607** | **1.83×** | 0.85 | 0.79 |
| L3 greedy | +1.009 | 2.74× | 2.73 | 0.73 |
| L3 NJ | +0.643 | 1.90× | 1.32 | 0.77 |
| L3 UPGMA | +0.708 | 2.03× | 1.47 | 0.78 |
| L3 ensemble | +0.642 | 1.90× | 1.51 | 0.67 |

**Branch-length miscalibration alone accounts for a 1.83-fold overestimate of the switching rate.** For NJ,
UPGMA and the ensemble, topology error adds almost nothing on top (0.04–0.10 in log units); only the greedy
reconstruction adds materially (+0.40).

This reverses the natural expectation. The field's effort goes into recovering topology, but in this benchmark
**the clock assumption is the larger problem.** Treating one resolved branch as one unit of time, when branches
in fact span variable numbers of divisions, biases the rate upward by about a factor of two before any
topological error is incurred.

What reconstruction *does* damage is precision: the IQR of the paired error rises from 0.85 at L2 to 1.3–2.7 at
L3. Reconstruction roughly doubles to triples the spread of the error while barely shifting its centre.

## Result 2 — star-likeness is biological, not a reconstruction artefact

This was the central question. Fraction of leaf pairs coalescing in the top 30% of the tree
("close-relative pairs", the pairs that carry switching-rate information):

| Recorder / genealogy | **true sampled tree** | greedy | NJ | UPGMA |
|---|---|---|---|---|
| cassiopeia_like, balanced | **0.016** | 0.012 | 0.019 | 0.014 |
| cassiopeia_like, bottlenecked | **0.019** | 0.012 | 0.019 | 0.012 |
| gestalt_like, balanced | **0.017** | 0.016 | 0.032 | 0.019 |
| gestalt_like, bottlenecked | **0.016** | 0.018 | 0.030 | 0.019 |

**The true sampled genealogy already has only 1.6–1.9% close-relative pairs**, and reconstructions give
1.2–3.2% — the same order, sometimes more. Normalised mean MRCA depth tells the same story (truth 0.24–0.26;
reconstructions 0.20–0.36).

So the scarcity of informative pairs that Milestone 4 found in the published trees is **a property of sampling
a modest number of cells from a large expanding population**, not an artefact of tree inference. Sampled cells
are mostly distant relatives because the population is large; no reconstruction method can change that.

**This is the answer to the Milestone 5 central question, and it is a negative one for algorithmic fixes.** The
missing information was never in the sample. Milestone 4's 0.70-decade penalty relative to the idealised
capacity law stands, and is attributable to sampling rather than to inference.

## Result 3 — topology recovery, and what capacity buys

Fraction of true clades recovered:

| Recorder | greedy | NJ | UPGMA |
|---|---|---|---|
| cassiopeia_like (25 sites, 72% saturated) | 0.03–0.04 | **0.56–0.60** | **0.57–0.60** |
| gestalt_like (9 sites, **97% saturated**) | 0.05–0.06 | 0.19–0.25 | 0.19–0.26 |

Higher recorder capacity roughly **doubles** topology recovery (0.56–0.60 against 0.19–0.26). The low-capacity
setting saturates at 97% within ten generations, so the recent divisions leave no mark and cannot be recovered
by any algorithm — Proposition V3 in practice.

The greedy solver performs poorly here, but **this should not be read as a verdict on Cassiopeia**: our greedy
is a simplified reimplementation without the priors and refinement steps the published solver uses. It is
included because it supplies reconstruction diversity for the ensemble, not as a benchmark of the authors'
software.

## Result 4 — the ensemble helps precision, not bias
Averaging the likelihood over three reconstructions gives the same median bias as the best single tree (+0.642
against +0.643 for NJ) but reduces the fraction of replicates that overestimate (0.67 against 0.73–0.78). In
the real data the effect is larger: the ensemble halves the apparent rate-heterogeneity estimate
(`milestone5_tree_reconstruction.md`). Ensembling therefore removes reconstruction *noise* but cannot remove
the clock bias, which is common to every tree.

## Result 5 — the information-processing chain is respected
Proposition V1 requires information to be non-increasing along
$(G_1,X_S) \to (B,X_S) \to (\hat G,X_S)$. Measured as inverse variance of $\log\hat k$, the observed ratios
against L1 are 1.24 (L2), 1.16–1.24 (single reconstructions) and 1.62 (ensemble) — all ≥ 1, i.e. no level shows
*more* information than the true tree. The apparent "gain" for the ensemble is a variance reduction accompanied
by unchanged bias, which Proposition V2 warns must not be read as recovered information.

## Consequences for practice
1. **Calibrate branch lengths, or model them.** Converting resolved branches to time is the largest single
   error. A recorder with a known edit rate, or an explicit model of divisions per resolved branch, would
   remove more bias than any improvement in tree search.
2. **Do not expect better reconstruction to rescue switching rates.** The informative pairs are absent from the
   sample itself.
3. **Sample densely within clades, not broadly across the tumour**, if switching kinetics is the target — the
   opposite of what maximises occupancy precision.
4. **Use an ensemble** where multiple reconstructions exist; it costs nothing and removes reconstruction noise.

## Limits
Synchronous divisions; independent per-site editing; three reconstruction algorithms, two of them
distance-based; a simplified greedy standing in for a published solver. The estimator itself is biased upward
at low rates even at L1 (see `milestone5_uncertainty_calibration.md`), so the paired design is what carries the
conclusions, not the absolute values.
