# Data provenance and redistribution policy

## Summary

This repository redistributes **no primary data and no third-party source code**. It
contains original analysis code, derived summary tables computed by that code, and
documentation. Everything else is retrieved by the user from its original source.

## Primary datasets

### Simeonov et al. 2021 — pancreatic adenocarcinoma, macsGESTALT

- **GEO: GSE173958** — 10x matrices.
- **Mendeley: doi:10.17632/t98pjcd7t6.1** — the authors' processed deposit, which supplies
  `colData`, clone assignments and the clone edgelists used to build trees.
- Samples used: mouse 1 GSM5283482/83/84/85/86/87, mouse 2 GSM5283488/89/90/91.

### Quinn et al. 2021 — lung adenocarcinoma xenograft, Cassiopeia

- **GEO: GSE161363** — 10x matrices and the authors' published reconstructions.
- Four animals (5k, 30k, 10k, 100k). Three reconstructions per lineage group were deposited
  for the 5k animal; the main analysis uses `hybrid_priors` only, except where the
  three-reconstruction ensemble is explicitly invoked.
- Barcode identifiers in the 10k animal carry a site prefix absent from the others; the
  loader selects whichever convention matches more cells.

## Derivation chain

```
GSE173958 10x matrices
  └─ src/biological_validation/extract_panel.py   (streams, keeps a fixed ~90-gene panel)
       └─ panel_M{1,2}.csv.gz
            └─ experiments/simeonov_reanalysis/run_real.py
                 └─ m2_real_M{1,2}_cells.csv.gz        [WITHHELD, see below]
                      └─ all m2_/m3_ summary tables     [committed]

GSE161363 matrices + published trees
  └─ src/biological_validation/quinn_states.py    (cell-cycle score, fixed external gene set)
       └─ quinn_cells_{5k,30k,10k,100k}.csv.gz
            └─ experiments/quinn_analysis/run_quinn.py
                 └─ all m4_ summary tables             [committed]
```

A note on the Quinn state definition: the first principal component of expression was
**rejected** as a state because it correlates with mitochondrial read fraction (−0.67 to
+0.60 across animals). A cell-cycle score from a fixed external gene set is used instead.
Results under the PC1 definition are retained as `m4_quinn_*_pc1.csv` so the choice can be
audited.

## What is withheld, and why

| File | Size | Reason |
|---|---|---|
| `m2_real_M1_cells.csv.gz` | 4.8 MB | per-cell normalised expression, 10,608 cells × ~90 genes |
| `m2_real_M2_cells.csv.gz` | 5.7 MB | per-cell normalised expression, 13,834 cells × ~90 genes |

These are derived from another group's GEO submission. Although they contain only ~90 of
roughly 20,000 genes, they are still processed per-cell expression, so they are not
redistributed. `scripts/` regenerates them and
`data/manifests/withheld_derived_sha256.txt` lets you verify the result is byte-identical.

Also not redistributed: the Quinn published trees and the Simeonov clone edgelists. Both
belong to the original authors and are available at the accessions above.

## Third-party software

The benchmarked methods are **not vendored**. Install from source at these commits:

| Method | Repository | Commit | Declared licence |
|---|---|---|---|
| PATH v1.0 | github.com/landau-lab/PATH | `caffc5105678224f108b3dd61fea545f96235ffc` | MIT (LICENSE file present) |
| SCOUT | github.com/hrstuart/SCOUT | `c72befe0d052f3960c448a64bff776665ab8ea4c` | DESCRIPTION says "MIT + file LICENSE"; **no LICENSE file in the repository** |
| PASTRI (package) | github.com/yuxiaochen11/PASTRI | `67aefd96240f9d3a1900c6906d813962233ef38f` | DESCRIPTION says "MIT license"; **no LICENSE file in the repository** |
| PASTRI (paper scripts) | github.com/Yangwj01/PASTRI | `5ecc99dfc5b7ede11bd20bf45ae47a1883267d5e` | none stated |
| OUwie 3.0.3 | CRAN | — | GPL-2 \| GPL-3 |

Two of the three declare MIT in `DESCRIPTION` but ship no licence text. Because the terms
are not unambiguous, **none of their source is redistributed here**; only our own adapters
are, and those are MIT. Anyone reproducing the benchmark installs each method directly from
its authors under whatever terms they offer.

Note that OUwie is GPL. We call it as an ordinary R dependency, installed by the user; no
OUwie code is included in this repository, so this repository's MIT licence is unaffected.

## What this repository does contain

- Original analysis code (`src/`, `experiments/`, `benchmarks/src/`, `scripts/`) — MIT.
- Derived **summary** tables in `results/tables/` and `benchmarks/results/tables/`,
  computed by that code: variance components, test statistics, calibration rates,
  per-clone fits. No per-cell expression.
- Figures generated from those tables.
- Protocols, decision records and audits.
