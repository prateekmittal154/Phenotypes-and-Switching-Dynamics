# Data

**No primary data is redistributed in this repository.** Both datasets are public and are
retrieved by the user. See `docs/data_provenance.md` for the full chain.

## Accessions

| Dataset | Accession | What to download |
|---|---|---|
| Simeonov *et al.* 2021 (pancreatic, macsGESTALT) | **GSE173958** | 10x matrices for the samples listed below |
| Simeonov processed deposit | Mendeley **doi:10.17632/t98pjcd7t6.1** | `colData`, clone assignments, clone edgelists |
| Quinn *et al.* 2021 (lung xenograft, Cassiopeia) | **GSE161363** | 10x matrices and the authors' published trees |

Simeonov GEO samples used: mouse 1 — GSM5283482 (PT), GSM5283483 (SS), GSM5283484 (Met),
GSM5283485 (Liver), GSM5283486 (Lung), GSM5283487 (CTC); mouse 2 — GSM5283488 (PT),
GSM5283489 (Met), GSM5283490 (Liver), GSM5283491 (Lung).

## Expected layout after retrieval

```
data/
├── raw/            <- what you download (git-ignored)
├── processed/      <- output of the extract/derive scripts (git-ignored)
└── manifests/      <- checksums, committed
```

## Retrieval and preprocessing

```bash
# Simeonov: stream the 10x matrices, keep the fixed ~90-gene panel
python -I src/biological_validation/extract_panel.py data/raw/simeonov \
       data/raw/simeonov_mendeley data/processed/

# Quinn: derive the cell-cycle state from a fixed external gene set
python -I src/biological_validation/quinn_states.py data/raw/quinn data/processed/
```

`extract_panel.py` streams rather than loading whole matrices: the Simeonov matrices are
about 1.6 GB gzipped.

## Withheld derived files

Two tables produced by the Simeonov chain carry per-cell normalised expression for ~90
genes across 24,442 cells:

- `m2_real_M1_cells.csv.gz` (10,608 cells × 107 columns)
- `m2_real_M2_cells.csv.gz` (13,834 cells × 107 columns)

They are **not committed**, to avoid redistributing processed per-cell expression derived
from another group's GEO submission. The scripts above regenerate them. Verify with:

```bash
sha256sum -c data/manifests/withheld_derived_sha256.txt
```

Every other table under `results/tables/` is a summary statistic computed by this code and
is committed.

## Manifests

`manifests/` holds the SHA-256 manifests from the primary study and the benchmark, carried
over so that a re-run can be checked against the frozen originals.
