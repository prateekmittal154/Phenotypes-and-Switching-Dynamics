# cancer-lineage-state-dynamics

Computational analysis of phenotypic inheritance, state occupancy and switching dynamics in
cancer lineage tracing.

---

## 1. What this is, and what it does

Single-cell lineage recorders promise to separate two explanations for a tumour changing
phenotype: **selection** among clones that already differ, and **induction** of new states
within clones. Both are usually summarised by one measurement — how alike related cells
are. This code shows that measurement is the sum of two different quantities, and that they
are not equally recoverable.

For a two-state process on a lineage, the within-clone covariance of states between two
cells at ancestral distance `L` decomposes exactly as

```
Cov(X_i, X_j | L)  =  E_c[ pi_c (1 - pi_c) exp(-2 k_c L) ]  +  Var_c(pi_c)
                      \_______ state memory _______/          \_ occupancy _/
```

- **`V_pi = Var_c(pi_c)`** — across-clone variance in how much of the time a state is held.
  It is the plateau the covariance curve approaches, and it is **invariant to the switching
  rate**.
- **`V_log k`** — across-clone variance in how fast cells interconvert. All the information
  about it sits in the decaying term, which is informative only for cell pairs with
  `k·L` near 0.43.

The software provides: the estimators for both quantities, the identifiability window that
says which rates a given tree can resolve, a birth–death–switching simulator with a full
observation chain (true genealogy → sampled → barcoded → reconstructed), calibration
machinery, and adapters that run **PATH**, **PASTRI** and the **OUwie** inference engine
on identical inputs for comparison.

**Headline result.** Clone-associated occupancy is recoverable and agrees across four
independent methods (median rank agreement 0.85 on real data). Switching rates are not:
on identical trees the three rate estimators disagree by more than an order of magnitude
and bracket the truth from both sides.

## 2. Installation

Python 3.12, plus R 4.4 only if you want to run the benchmarks.

```bash
git clone https://github.com/<owner>/cancer-lineage-state-dynamics.git
cd cancer-lineage-state-dynamics

conda env create -f environment.yml
conda activate cldyn

pip install -e .          # editable; puts the repo root on sys.path
```

The analysis modules import as `src.*` and `experiments.*` and expect the **repository root
as the working directory**. This is deliberate: the code is the archived, verified version
used for the paper, and renaming its import paths would have meant editing ~50 files that
were validated as they stand. `pip install -e .` makes those imports resolve from anywhere;
running from the repo root also works without installing.

Verify:

```bash
python -c "import src.lineage_inference.ctmc as m; print(m.__file__)"
pytest -q
```

R dependencies for the benchmarks are in `benchmarks/environments/README.md`.

## 3. Minimal working example

Simulate a 12-clone lineage where clones differ only in occupancy, then recover it:

```bash
python -m scripts.minimal_example
```

Or in Python:

```python
import numpy as np
from src.branching.recorder import simulate_genealogy, induced_subtree
from src.lineage_inference.ctmc import VarTree, fit_var

rng = np.random.default_rng(0)
par, gen, live = simulate_genealogy(10, rng)          # synchronous genealogy
take = rng.choice(live, 48, replace=False)            # sample 48 cells
sub_par, sub_elen, leafmap, _ = induced_subtree(par, take)

# two-state CTMC with occupancy 0.7, switching rate 0.3 per generation
E = np.exp(-0.3); st = np.empty(len(par), int); st[0] = 1
for i in range(1, len(par)):
    p = st[par[i]]
    st[i] = int(rng.random() < (0.7 * (1 - E) if p == 0 else 0.7 + 0.3 * E))

leaf = np.full(len(sub_par), -1)
for new_i, orig in leafmap.items():
    leaf[new_i] = st[orig]

print(fit_var(VarTree(sub_par, leaf, elen=sub_elen)))
# -> {'pi': ~0.7, 'k': ~0.3, 'loglik': ..., 'k_at_bound': False}
```

`k_at_bound` is the flag to watch: when it is `True` the tree does not resolve the rate and
the estimate should not be read as a measurement.

## 4. Datasets: retrieval and preprocessing

**No primary data is redistributed here.** Both datasets are public:

| Dataset | Accession | Contents |
|---|---|---|
| Simeonov *et al.* 2021, pancreatic, macsGESTALT | **GSE173958** + Mendeley `doi:10.17632/t98pjcd7t6.1` | 10x matrices, clone edgelists |
| Quinn *et al.* 2021, lung xenograft, Cassiopeia | **GSE161363** | 10x matrices, published trees |

Retrieval and the exact preprocessing chain are in **`docs/data_provenance.md`**. In brief:

```bash
# Simeonov: stream the 10x matrices and keep the fixed gene panel
python -I src/biological_validation/extract_panel.py <raw_dir> <mendeley_dir> data/processed/

# Quinn: derive the cell-cycle state from a fixed external gene set
python -I src/biological_validation/quinn_states.py <quinn_raw_dir> data/processed/
```

Two derived tables (`m2_real_M1_cells.csv.gz`, `m2_real_M2_cells.csv.gz`) carry per-cell
normalised expression for ~90 genes and are **withheld** rather than redistributed. The
scripts above regenerate them; `data/manifests/withheld_derived_sha256.txt` holds their
SHA-256 so you can confirm you reproduced them byte-for-byte.

## 5. Reproducing the mathematical experiments

All synthetic results, from the repository root:

```bash
bash scripts/run_all_studyd.sh
```

Or individually:

```bash
python -I experiments/run_all.py                                   # bounds, scenarios, design
python -I experiments/lineage_identifiability/run_e8.py            # covariance signatures
python -I experiments/lineage_identifiability/run_e9.py            # mechanism recovery
python -I experiments/lineage_identifiability/run_e10_focus.py     # occupancy vs rate
python -I experiments/sampling_design/fisher_audit.py              # kL* = 0.43
python -I experiments/sampling_design/window_recompute.py          # identifiability window
python -I experiments/sampling_design/full_tree_information.py --mc  # pair-sum inflation
python -I experiments/calibration/run_e14.py                       # interval calibration (~1 h, 12 cores)
python -I experiments/sampling_design/run_e15b.py                  # recorder capacity
python -I experiments/sampling_design/capacity_verification.py results/tables
python -I experiments/reconstruction/run_e17.py                    # observation hierarchy
python -I experiments/reconstruction/run_e20.py                    # reconstruction-aware intervals
python -I experiments/sampling_design/run_e21_sampling.py          # sampling design
```

Outputs land in `results/tables/`. Seeds are fixed inside each script.

## 6. Reproducing the biological analyses

These need the data retrieved in §4 and are **not** in `run_all_studyd.sh`, because the
analysis server cannot reach NCBI.

```bash
python -I experiments/simeonov_reanalysis/run_real.py        # pre-registered analysis
python -I experiments/simeonov_reanalysis/run_controls.py    # depth / technical controls
python -I experiments/simeonov_reanalysis/run_models.py      # nested CTMC comparison
python -I experiments/simeonov_reanalysis/run_m3.py          # per-clone estimation
python -I experiments/quinn_analysis/run_quinn.py            # Quinn occupancy and rates
python -I experiments/quinn_analysis/quinn_controls.py       # technical control
python -I experiments/quinn_analysis/quinn_primary_only.py   # depth-adjusted analysis
python -I experiments/reconstruction/run_e22_realboot.py     # dataset-specific bootstrap null
```

## 7. Running the benchmark comparisons

PATH, PASTRI and the OUwie engine are **not vendored**. Install them from their own
repositories at the exact commits used here:

```r
remotes::install_github("landau-lab/PATH",      ref = "caffc5105678224f108b3dd61fea545f96235ffc")
remotes::install_github("yuxiaochen11/PASTRI",  ref = "67aefd96240f9d3a1900c6906d813962233ef38f")
install.packages("OUwie")   # 3.0.3 used here
```

Then:

```bash
# verify each method reproduces its authors' own example before comparing
Rscript benchmarks/src/adapters/reproduce_examples.R $PWD/benchmarks

# generate benchmark data and run all four methods
python -I benchmarks/src/simulations/scenarios_bench.py . benchmarks/data/sim_ground S3 25 7000
bash benchmarks/src/adapters/run_gt_all.sh      # ground truth + branch-time calibration
bash benchmarks/src/adapters/run_phase2.sh      # tree geometry + real data
bash benchmarks/src/adapters/run_phase3.sh      # reconstruction uncertainty
bash benchmarks/src/adapters/run_phase4.sh      # held-out prediction

# score and draw
bash benchmarks/src/finalise.sh
```

**On SCOUT:** only its **OUwie inference engine** was executed, through
`benchmarks/src/adapters/scout_fn.R`. The full `SCOUT()` pipeline — Pagel's-lambda
screening, the EM step, regime assignment and hybrid model selection — was **not run**, and
no claim here describes SCOUT's complete pipeline. See
`benchmarks/reports/method_reproduction.md`.

One documented substitution: PASTRI imports `ggtree` for a single `fortify()` call, which
could not be built in our environment; it is supplied from `tidytree` in
`benchmarks/src/adapters/pastri_load.R`. The estimator itself is the authors' unmodified
code. See `benchmarks/reports/deviations.md`.

## 8. Regenerating the figures

```bash
export PYTHONPATH=<path to a site dir providing cycler>   # only if matplotlib deps are split
python experiments/figures.py
python experiments/figures_m2.py results/tables results/figures
python experiments/figures_m3.py results/tables results/figures
python experiments/figures_m4.py results/tables results/figures
python experiments/figures_m5.py results/tables results/figures

python -I benchmarks/src/metrics/make_bench_figures.py $PWD/benchmarks
```

Every panel reads a CSV in `results/tables/` or `benchmarks/results/tables/`; no figure
contains a hand-entered number. Panel-to-table provenance for the benchmark figures is in
`benchmarks/reports/`.

## 9. Hardware and expected cost

Developed on a 24-core Linux server, 30 GB RAM. Nothing needs a GPU.

| Step | Cost |
|---|---|
| `pytest -q` | < 1 min |
| `scripts/minimal_example.py` | < 10 s |
| Milestone 1–3 synthetic experiments | ~10 min |
| `run_e14.py` interval calibration | ~1 h on 12 cores |
| `run_e17.py` / `run_e20.py` reconstruction | ~1–2 h each |
| Full `run_all_studyd.sh` | ~4 h |
| Benchmark ground truth (4 methods × 200 replicates × 2 conventions) | ~4 h on 10 cores |
| Benchmark real data (116 clone trees) | ~1 h |
| Peak memory | < 4 GB per process |

The benchmark adapters are single-threaded by design; set
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` so that *N* parallel jobs really use *N* cores.

## 10. Limitations and provenance

**What the analysis supports.** Clone-associated occupancy heterogeneity is measurable,
robust to the reconstruction used, and reproduced by four independent methods.

**What it does not.** Per-clone switching rates. In every published dataset examined the
fitted rates lie outside the band the trees resolve, and the four methods disagree. The
estimator here is not more accurate than the alternatives on rates — it **over-disperses**
where PATH and PASTRI compress. Its contribution is the diagnostic: an interval, an
explicit identifiability window, and a boundary flag.

Also not claimed: any mechanism for why clones differ; human relevance; that better
reconstruction would recover the rates; any gene-level or spatial validation.

Full detail:

- `docs/methodology.md` — the model, estimators and identifiability results
- `docs/reproducibility.md` — environment, seeds, checksums, what reproduces and what does not
- `docs/data_provenance.md` — accessions, retrieval, derivation chain, redistribution policy
- `protocols/deviations.md` — every logged deviation in the primary study
- `benchmarks/reports/deviations.md` — every logged deviation in the benchmark, including
  two bugs found and fixed in our own code
- `reports/` — the milestone decision records and claim-to-evidence audits

## Citation

See `CITATION.cff`. The manuscript is in preparation; this repository is the code release.

## Licence

MIT — see `LICENSE`. Third-party methods are not redistributed; install them from their own
repositories under their own licences.
