# Deviations and judgement calls

Everything here changed what was run or how it was interpreted. Nothing is omitted because
it was inconvenient.

---

**B1. ggtree could not be installed; PASTRI's `fortify()` is supplied from tidytree.**
The Bioconductor ggtree build for R 4.4 calls `ggplot2:::check_linewidth` (removed in
ggplot2 4.x) and `ggplot2::is_ggplot` (added in 4.x), so it fails against both 3.5.1 and
4.0.4; the GitHub build needs `ggiraph` -> `gdtools`, whose system libraries are absent.
PASTRI uses ggtree for one call, `fortify(tree)`, and reads only the `parent` and `node`
columns. We supply that from `tidytree` (same authors; ggtree's `fortify.phylo` is built on
it) and source PASTRI's three estimator files verbatim. The adjacency construction and the
Frobenius-norm optimisation are unmodified authors' code. See `src/adapters/pastri_load.R`.
*Effect:* none expected on results; the substituted function is a data-frame conversion.

**B2. ggplot2 pinned to 3.5.2.** Required because PATH imports `fgsea`, which needs
>= 3.5.2. Does not affect any computation used here.

**B3. SCOUT is represented by its OUwie engine, not the full `SCOUT()` pipeline.**
`SCOUT()` adds Pagel's-lambda screening, an EM step, regime assignment and a hybrid
selection rule, and its `NAMESPACE` imports `TedSim` (a GitHub simulation dependency).
We run `OUwie::OUwie` with `OU1`/`BM1` directly, per clone, which is the model-fitting step
SCOUT wraps. *Effect:* SCOUT's results here are a lower bound on what the full pipeline
would do. Any statement about SCOUT is scoped to the engine, and we do not claim its full
pipeline behaves the same way. `runSCOUT.identifiability()` was not run.

**B4. PASTRI's `fi_depth` is swept, not fixed.** The authors use cumulative windows of
0.4–0.8 in their own applications and state that including more distal adjacencies
overestimates rates. Fixing a single value would have been an arbitrary choice that could
flatter or penalise the method. We sweep {0.2, 0.4, 0.6, 0.8, 1.0}, report the whole sweep
(Figure B3A) and use 0.4 for headline numbers, matching the authors' C. elegans setting.
*Effect:* this is more favourable to PASTRI than a naive `fi_depth = 1.0`, which our pilot
showed drives the implied switch probability to its boundary.

**B5. PASTRI returns probabilities, not rates.** Its `T` is a per-cell-cycle transition
probability matrix with rows summing to 1. For a two-state chain the off-diagonal sum over
one generation is `1 - exp(-k)`, so we convert with `k = -log(1 - p)`. Cases where `p`
reaches 1 give an infinite rate; these are counted as "no estimate" rather than dropped
silently. *Effect:* without this conversion PASTRI and the rate-based methods would have
been compared on different scales.

**B6. Clone size was raised from 24 to 48 cells after a pilot.** The first pilot used 12
clones x 24 cells. PASTRI is designed for trees with hundreds to thousands of tips, and at
24 cells per clone its depth resolution was visibly poor. We regenerated all simulations at
48 cells per clone (576 per replicate) before running anything that is reported.
*Effect:* more favourable to PASTRI and to every method. Benchmark II varies depth and
sampling further so the size dependence is visible rather than hidden.

**B7. A monomorphic clone breaks PATH.** When every cell of a clone is in one state,
`catMat()` returns a single column and `PATH_inf()` cannot form a 2x2 matrix; the first run
of scenario S1 crashed on this. We added a guard that records "monomorphic clone" and
continues. *Effect:* this is why PATH's no-estimate rate in S1 is 0.130 rather than a crash.
It is a property of the method on near-frozen data, and is reported as such.

**B8. The two adapters labelled `scenario` differently and wrote `true` with different
quoting.** The R adapter used the replicate directory name and quoted strings; the Python
adapter used the scenario prefix and wrote `true` unquoted, which pandas read back as a
boolean. The first merge therefore matched nothing. Fixed by joining on the replicate
directory and normalising the variant column in all three analysis scripts. *Effect:* a
bug in the analysis layer only; caught before any result was recorded, and no reported
number predates the fix.

**B9. The method-reproduction script was trimmed from 60 to 40 genes.** The first run fitted
60 genes x 2 models with multithreaded BLAS and consumed about eight cores for over an hour,
starving the benchmark jobs. It was stopped and rerun single-threaded with 20 genes per
generating class. *Effect:* the model-selection accuracy (0.625) rests on 40 genes rather
than 60. The direction of the result — all OU genes recovered, most BM genes called OU — is
unambiguous at this count.

**B10. Study D is the method under test and its authors wrote this benchmark.** Mitigations
are listed in `protocols/evaluation_rules.md`. The most important are that every competitor
was run from its authors' published code and first reproduced on the authors' own example,
and that two scenarios (S2 rate-only, S5 time-varying induction) generate data outside
Study D's own time-homogeneous likelihood family. Results unfavourable to Study D — notably
that it produces the **largest** spurious rate heterogeneity of the three binary methods —
are reported in the main results, not buried.

**B11. A bug in our own real-tree exporter collapsed four Simeonov M1 clones to a single
tip, and was caught and fixed before any real-data result was recorded.**
`prune_to_labelled()` in `src/simulations/export_real_trees.py` selected the tree root as
`min(needed)`. That idiom is correct in the Study D simulator, where node indices are
assigned in birth order so the lowest index is the root, but the Simeonov nodes are indexed
by **sorted name**, so the lowest index is the alphabetically first node. For four of the
eight M1 clones the alphabetically first node was not the root, and the pruned tree
collapsed to a single leaf: `tree_unit.nwk` was 28 bytes, e.g.
`Blood_GTCATTTCATTCACAG:1.0;`, for clones carrying 45, 72, 105 and 238 cells.

Detected because `ape::read.tree` refused the degenerate file and the Simeonov M1
competitor job halted; Study D's own parser had silently skipped those clones under its
`len(cells) < 6` guard, so the failure would otherwise have shown up only as four missing
clones.

Fixed by selecting the root as the node with `parent < 0`, and by adding two guards that
raise rather than write: `newick_from_parent` refuses a forest, and `write_one` refuses a
tree with fewer than two labelled tips. All 116 clone trees were re-exported and every
real-data run was repeated from scratch. The faulty export is retained at
`data/real_trees_BUGGED_degenerate/` for audit. *Effect:* no reported real-data number
comes from the faulty trees; the four recovered M1 clones carry 615 of the animal's cells
between them.

**B12. SCOUT failed on every real clone tree at first pass, because of unary nodes.**
The published Simeonov and Quinn trees contain single-child nodes, including at the root.
`ape::multi2di` resolves polytomies but does not remove unary nodes, so `OUwie` rejected
every real tree with *"'phy' is not rooted and fully dichotomous"* and the first real-data
pass recorded `SCOUT_alpha = NA` for all 116 clones. Fixed by calling
`ape::collapse.singles()` before `multi2di()` in `run_ouwie_clone`, with an explicit
`is.rooted`/`is.binary` guard afterwards and a re-alignment of the trait vector to the
possibly reordered tips. The OUwie fit was then re-run over all 116 real clone trees
(`src/adapters/run_scout_real.R`) and merged into the real-data analysis; PATH, PASTRI and
Study D were not re-run, since nothing in their path changed. To keep the two code paths
from drifting, the function is now defined once in `src/adapters/scout_fn.R` and sourced
by both scripts. *Effect:* without this, SCOUT would have been absent from the real-data
comparison entirely. The simulated trees were unaffected (SCOUT returned a finite alpha for
600/600 held-out clones and 4,800/4,800 ground-truth clones).

---

## Not done

**N1. Group D held-out prediction was run on two scenarios, not all eight.** S3 and S4,
50 replicates. Enough to compare predictive skill across all four methods; not enough to
claim the ordering holds under every mechanism.

**N2. Benchmark IV uses one recorder parameterisation.** Cassiopeia-like only
(30 sites, edit rate 0.06, 40 alleles, 10% dropout, 10% missing). A GESTALT-like
parameterisation is implemented in `src/simulations/reconstruction_bench.py` but was not
run for want of time.

**N3. No method was given tuned hyperparameters.** Each was run at its documented defaults
and published settings. It is possible that tuning would improve any of them, and we make
no claim about their best achievable performance.

**N4. PASTRI's permutation test and SCOUT's identifiability utilities were not run.**
Both exist and both would add to the uncertainty picture. Neither is required for the
estimand comparisons reported here.
