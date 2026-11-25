# Reproducibility

## Environment

The analysis ran on Linux, Python 3.12 with numpy 2.5, scipy 1.18, pandas 3.0. The pinned
environment is `environment.yml`; `pyproject.toml` declares looser floors so the code also
installs on current releases.

```bash
conda env create -f environment.yml && conda activate cldyn
pip install -e .
pytest -q
```

The benchmark layer additionally needs R 4.4 with PATH, PASTRI and OUwie; see
`benchmarks/environments/README.md`. Two environment facts cost real time to discover and
are worth repeating:

- The conda compilers must be on `PATH` or every R source build fails with
  `x86_64-conda-linux-gnu-c++: not found`.
- `fgsea`, which PATH imports, needs ggplot2 >= 3.5.2. Pinning 3.5.1 to satisfy `ggtree`
  breaks PATH.

## Seeds and determinism

Every experiment script sets its seed explicitly at the top and passes a
`numpy.random.default_rng` through the call chain; none relies on global state. The
simulation generators take the seed on the command line and record it in each replicate's
`truth.json`. The estimators are deterministic given a tree.

## What reproduces exactly, and what does not

**Exactly.** All synthetic experiments, from `scripts/run_all_studyd.sh`. The benchmark
suite, from the four phase drivers. Every figure, from the deposited tables.

**Requires data retrieval first.** The biological analyses. Both datasets are public but
must be downloaded; the analysis server used for this work cannot reach NCBI, so these
steps were run separately. See `docs/data_provenance.md`.

**Known to differ.** Nothing in the reported results depends on hardware, but note that
`results/tables/` holds the frozen outputs of the runs reported in the manuscript. A re-run
regenerates them; compare against `data/manifests/` rather than assuming byte equality
across BLAS versions.

## Checksums

`data/manifests/` carries the SHA-256 manifests from the primary study and the benchmark.
The Milestone 6 manifest, which covers the final state of the primary analysis, verified
644/644 files at the time of release. Two frozen protocol manifests
(`FREEZE_SHA256_m1.txt`, `FREEZE_SHA256_m2_simeonov.txt`) verify as well.

Earlier milestone manifests show a handful of expected mismatches — `deviations.md`,
`ctmc.py`, `perclone.py`, `recorder.py` — all files legitimately revised in later
milestones. The Milestone 6 manifest supersedes them.

## Cost

| Step | Wall clock |
|---|---|
| `pytest -q` | < 1 min |
| minimal example | < 10 s |
| synthetic experiments, Milestones 1–3 | ~10 min |
| `run_e14.py` calibration | ~1 h, 12 cores |
| `run_e17.py`, `run_e20.py` | ~1–2 h each |
| full `run_all_studyd.sh` | ~4 h |
| benchmark ground truth | ~4 h, 10 cores |
| benchmark real data | ~1 h |

Set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` before running parallel jobs; otherwise
each R process grabs every core and *N* jobs run slower than one.

## Corrections made during the work

Both studies log every deviation. Two classes are worth flagging to anyone rebuilding this:

**In the primary study** (`protocols/deviations.md`, 29 entries) — notably a
Robinson–Foulds comparison that reported zero clade recovery because sampled and
reconstructed leaves were indexed in different label spaces, and a naive global
branch-scale correction that was implemented, found to *worsen* calibration, and retained
as a reported negative result rather than removed.

**In the benchmark** (`benchmarks/reports/deviations.md`, 12 entries) — including two bugs
in our own code, both found and fixed before any result was recorded:

- B11: the real-tree exporter selected the root as `min(needed)`, correct in the simulator
  where nodes are indexed in birth order but wrong for Simeonov nodes indexed by sorted
  name. Four clone trees collapsed to a single tip. All 116 trees were re-exported and
  every real-data run repeated.
- B12: `ape::multi2di` does not remove unary nodes, so OUwie rejected all 116 real trees.
  Fixed with `collapse.singles()`; the OUwie call is now defined once in
  `benchmarks/src/adapters/scout_fn.R` and shared, so the two code paths cannot drift.

## Scope of the claims

Stated so results are not read more broadly than they support.

1. **Switching-rate intervals are regime-dependent.** Pooled coverage mixes an identifiable
   and a non-identifiable regime. Coverage is reported separately inside and outside the
   window, and for real data only the window each dataset resolves is reported, with no
   point estimate below it.
2. **The benchmark ran the original PATH and PASTRI implementations and the OUwie engine.**
   It did **not** run the full `SCOUT()` pipeline (lambda screening, EM, regime assignment,
   hybrid selection), PASTRI's permutation test, or a GESTALT-like recorder in the
   reconstruction benchmark.
3. **The sampling-design result is split.** The informative-pair-fraction gain is robust at
   40 observations per strategy; the coverage ordering has overlapping Wilson intervals and
   is reported as suggestive.
4. **Branch times are internal.** Ultrametric times come from each tree's own depth
   structure. No external molecular clock was imposed.
5. **Not claimed.** Any mechanism for why clones differ; human relevance; that better
   reconstruction would recover the rates; any gene-level or spatial validation. Nine
   held-out gene programmes produced no association with clone occupancy that replicated
   across animals, with three changing sign.

## Known unresolved cross-references

`python tools_link_check` style validation over the 68 markdown files leaves eleven
references that do not resolve. All are expected; none affects running the code.

| Reference | Why it does not resolve |
|---|---|
| `theory/counterexamples/C1.md`, `C7.md`, `../counterexamples/C2.md` | Pre-existing drift in the archived theory notes: the counterexamples were consolidated into `C1_to_C5.md` and `C6_C8.md`. The prose is left as written rather than edited after the fact. |
| `results/tables/e10_occupancy_vs_rate.csv` | Same drift: the table shipped as `e10_focus.csv`. |
| `data/manifests/dataset_suitability_matrix.csv` | Referenced in an early feasibility report; the file was never produced. |
| `manuscript/journal_targets.md`, `submission/external_assets_needed.md` | Live in the manuscript and submission packages, which are deliberately **not** part of this code release. |
| `R/xcor.R` (PATH), `R/LCA_depths_calculation.R` (PASTRI) | Files inside the third-party repositories, which are not vendored here. Install the packages to read them; see `docs/data_provenance.md`. |
| `data/real_trees_BUGGED_degenerate/` | Generated data retained on the analysis machine as audit evidence for deviation B11; not shipped because it is regenerable and large. |

The first four are inherited from the canonical research directories and were left intact
so the archived records read as they did when written.
