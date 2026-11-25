"""Build and execute notebooks/milestone2_results.ipynb."""
import os, nbformat as nbf
from nbclient import NotebookClient
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = [
    md("# Study D, Milestone 2\nDecision: **GO**, scoped to a methods-and-re-analysis paper. See "
       "`reports/milestone2_decision.md`."),
    code("import json, sys, pandas as pd, numpy as np\nsys.path.insert(0, '..')\npd.set_option('display.width', 220)"),
    md("## Synthetic ground truth: mechanism recovery (E9, primary tumour)\n"
       "`naive` is the variance of clone-level state fractions, the statistic a clone analysis reports by default."),
    code("d = pd.read_csv('../results/tables/e9_mechanism_recovery.csv')\n"
         "d[d.organ == 0].groupby('mech').agg(truth_B=('truth_B','median'), B=('B','median'),\n"
         "    naive=('clone_frac_var','median'), frac_CI_excludes_0=('B_lo', lambda s: (s>0).mean()),\n"
         "    reps=('B','size')).round(4)"),
    md("### Covariance curves by MRCA depth (the signature)"),
    code("c = pd.read_csv('../results/tables/e9_curves.csv')\n"
         "c[c.organ==0].pivot_table(index='depth', columns='mech', values='site_centred').round(4)"),
    md("### Stress tests: misclassification and clone-label loss attenuate B, never inflate it"),
    code("s = pd.read_csv('../results/tables/e9_stress.csv'); s = s[s.organ==0]\n"
         "s.groupby(['eps','clone_loss','n_sample']).B.median().round(4).unstack()"),
    md("## Real data: Simeonov GSE173958\n### Where the metastatic EMT enrichment comes from"),
    code("pd.read_csv('../results/tables/m2_real_c3_composition.csv').round(3)"),
    md("### Technical control: many unrelated programmes shift as much as EMT"),
    code("t = pd.read_csv('../results/tables/m2_real_c1_technical_control.csv')\n"
         "t.pivot_table(index='feature', columns=['mouse','site'], values='within_clone_shift_sd').round(2)"),
    md("### State heritability (clone-stratified permutation)"),
    code("pd.read_csv('../results/tables/m2_real_c2_h2_permutation.csv').round(4)"),
    md("### The strict test: heritable switching propensity"),
    code("pd.read_csv('../results/tables/m2_real_decomposition.csv')[\n"
         "  ['mouse','site','n_clones','n_cells','k','A','B','B_lo','B_hi']].round(4)"),
    md("### Model comparison and held-out site prediction"),
    code("pd.read_csv('../results/tables/m2_real_model_comparison.csv')[\n"
         "  ['mouse','model','npar','loglik','aic','delta_aic']].round(2)"),
    code("h = pd.read_csv('../results/tables/m2_real_heldout_site.csv')\n"
         "h.groupby('mouse')[['err_M0','err_M1']].mean().round(4)"),
    code("json.load(open('../results/metrics/m2_real_metrics.json'))"),
]
nb = nbf.v4.new_notebook(cells=cells)
NotebookClient(nb, timeout=900, kernel_name="python3",
               resources={"metadata": {"path": os.path.join(ROOT, "notebooks")}}).execute()
nbf.write(nb, os.path.join(ROOT, "notebooks/milestone2_results.ipynb"))
print("notebook ok")
