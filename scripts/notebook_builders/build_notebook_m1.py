"""Build and execute notebooks/milestone1_results.ipynb (needs nbformat/nbclient from Study A env/site)."""
import os, nbformat as nbf
from nbclient import NotebookClient

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = [
    md("# Study D, Milestone 1: results walkthrough\nDecision: **REFORMULATE** (see `reports/milestone1_decision.md`)."),
    code("import json, sys, numpy as np, pandas as pd\nsys.path.insert(0, '..')\n"
         "from src.bounds.core import tau_min, tau_min_lp, tv, support_mismatch, tau_lower_robust\n"
         "m = json.load(open('../results/metrics/milestone1_metrics.json'))"),
    md("## P2: closed form vs independent LP (E1)"),
    code("pd.Series(m['E1'])"),
    md("## C1: one (P, Q), many mechanisms"),
    code("P = np.array([.9, .09, .01]); Q = np.array([.4, .3, .3])\n"
         "pd.DataFrame({'Lambda': L, 'tau_min': tau_min(P, Q, L), 'LP': tau_min_lp(P, Q, L)} for L in [1, 2, 5, 10, 30, 67.5, 100])"),
    md("## Nontriviality (E2): frozen family and the biologically shaped family"),
    code("d = pd.read_csv('../results/tables/e2_nontriviality_perturbation_family.csv')\n"
         "print(m['E2'])\n"
         "d[d.lam_true <= d.lam].groupby('lam').agg(n=('i', 'size'), frac_positive=('tau_min', lambda x: (x > 1e-9).mean()), tau_true=('tau_true', 'median'))"),
    md("## Scenarios S1-S8 (E4)"),
    code("e = pd.read_csv('../results/tables/e4_scenarios.csv')\n"
         "e[['scenario', 'n', 'lam', 'tau_true', 'tau_min_pop', 'plug_pos', 'tv_pos', 'sm_pos', 'robust_pos', 'robust_cover', 'robust_mean']].round(3)"),
    code("pd.read_csv('../results/tables/e4_equivalence_twins.csv')"),
    code("pd.Series(m['E4_S5'])"),
    md("## Design (E5) and founder bottleneck (E6)"),
    code("pd.read_csv('../results/tables/e5_design_required_n.csv').pivot_table(index=['q', 'rho'], columns='lam', values='n_required_cells')"),
    code("pd.read_csv('../results/tables/e5_regional_bias.csv')"),
    code("pd.read_csv('../results/tables/e6_single_met_neutral.csv').pivot(index='nf', columns='lam', values='false_induction')"),
    code("pd.read_csv('../results/tables/e6_fluctuation_tests.csv').round(4)"),
    md("## Dynamics (E7)"),
    code("pd.read_csv('../results/tables/e7_dynamics_fits.csv').groupby(['design', 'ncell'])[['a', 'c', 'db']].median().round(3)"),
]
nb = nbf.v4.new_notebook(cells=cells)
NotebookClient(nb, timeout=600, kernel_name="python3", resources={"metadata": {"path": os.path.join(ROOT, "notebooks")}}).execute()
nbf.write(nb, os.path.join(ROOT, "notebooks/milestone1_results.ipynb"))
print("notebook ok")
