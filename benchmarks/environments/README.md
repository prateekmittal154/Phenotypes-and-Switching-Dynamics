# Environment

The conda R environment used for every run is **not included in the archive** (2.2 GB).
It is rebuilt from the scripts here, in this order:

```bash
B=<path to this package>
conda create -y -p $B/environments/renv -c conda-forge \
  r-base=4.4 r-ape r-devtools r-matrix r-rcpp r-rcpparmadillo r-castor r-phytools r-nlme r-mass

conda install -y -p $B/environments/renv -c conda-forge \
  r-dplyr r-data.table r-tidyverse r-reshape2 r-nloptr r-corpcor r-magrittr r-stringr \
  r-rlang r-future r-future.apply r-progressr r-foreach r-doparallel r-remotes r-argparse \
  r-patchwork r-tidyr r-plyr zlib libpng mpfr gmp

# The conda compilers must be on PATH or every source build fails with
#   x86_64-conda-linux-gnu-c++: not found
export PATH=$B/environments/renv/bin:$PATH

Rscript $B/environments/install_deps.R    # CRAN packages for PATH, SCOUT, PASTRI
Rscript $B/environments/install_deps3.R   # png, Rmpfr, geiger, paleotree, OUwie
Rscript $B/environments/install_bioc.R    # fgsea
```

Two pins matter.

- **ggplot2 3.5.2.** `fgsea`, which PATH imports, requires >= 3.5.2. An earlier attempt
  pinned 3.5.1 to satisfy ggtree and that broke fgsea, so PATH would not load.
- **ggtree is not installed and is not needed.** See `reports/deviations.md` B1: PASTRI's
  single `fortify()` call is supplied from `tidytree` in `src/adapters/pastri_load.R`, and
  the estimator itself is the authors' unmodified code.

Exact versions of everything used are in `results/tables/TableB2_software.csv`.

## Method source

`methods/` holds the four repositories as cloned, with `.git` removed to keep the archive
small. The commits actually used were:

| Repository | Commit |
|---|---|
| github.com/landau-lab/PATH (v1.0) | `caffc5105678224f108b3dd61fea545f96235ffc` |
| github.com/hrstuart/SCOUT | `c72befe0d052f3960c448a64bff776665ab8ea4c` |
| github.com/yuxiaochen11/PASTRI (package) | `67aefd96240f9d3a1900c6906d813962233ef38f` |
| github.com/Yangwj01/PASTRI (paper scripts) | `5ecc99dfc5b7ede11bd20bf45ae47a1883267d5e` |

## Python

Analysis and simulation use the Study D interpreter
(`Study_D_Metastatic_Adaptation_Identifiability/venv_shared/bin/python`): Python 3.12,
numpy 2.5, scipy 1.18, pandas 3.0. Figures additionally need the shared site directory at
`Study_A_Vascular_Ecosystem/env/site` for `cycler`, which `make_bench_figures.py` adds to
`sys.path` itself.
