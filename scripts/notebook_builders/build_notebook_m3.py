"""Build and execute notebooks/milestone3_results.ipynb."""
import os, nbformat as nbf
from nbclient import NotebookClient
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
T = "../results/tables/"
cells = [
    md("# Study D, Milestone 3\n"
       "Corrects the Milestone 2 interpretation, separates clone occupancy from switching rate, and derives the "
       "lineage sampling-design rule. See `reports/milestone3_decision.md`."),
    code("import json, os, sys, numpy as np, pandas as pd\nsys.path.insert(0, '..')\n"
         "pd.set_option('display.width', 220)\nT = '../results/tables/'\n"
         "have = lambda f: os.path.exists(T + f)"),
    md("## 1. The audit (U1): the plateau is invariant to the switching rate"),
    code("expit = lambda x: 1/(1+np.exp(-x)); logit = lambda p: np.log(p/(1-p))\n"
         "pi = np.array([0.2, 0.8])\n"
         "pd.DataFrame([{'k': k, 'plateau Var_c(pi)': pi.var(), 'A(0)': (pi*(1-pi)).mean()}\n"
         "              for k in (0.01, 1.0, 10.0)])"),
    md("### U3: a common induction shift rescales the plateau but cannot create one"),
    code("pd.DataFrame([{'logit shift': d, 'pi_c shifted': np.round(expit(logit(pi)+d), 3),\n"
         "               'Var (clones differ)': expit(logit(pi)+d).var(),\n"
         "               'Var (clones identical)': expit(logit(np.array([.5,.5]))+d).var()}\n"
         "              for d in (0, .5, 1, 2)])"),
    md("## 2. Separating occupancy from switching rate (E10-focus, 20 clones, 10 replicates)\n"
       "OCC and FROZE share an occupancy distribution and differ only in timescale. FROZE sits below the\n"
       "identifiability window, where the point estimates break but the likelihood-ratio tests stay calibrated."),
    code("f = 'e10_focus.csv' if have('e10_focus.csv') else 'e10_occupancy_vs_rate.csv'\n"
         "d = pd.read_csv(T+f)\n"
         "cols = [c for c in ['true_V_pi','V_pi','true_V_logk','V_logk','n_clones_k_identified',\n"
         "                    'lrt_pi_p','lrt_k_p'] if c in d.columns]\n"
         "display(d.groupby('mech')[cols].median().round(4))\n"
         "display(d.assign(rejects_shared_pi=d.lrt_pi_p<.05, rejects_shared_k=d.lrt_k_p<.05)\n"
         "         .groupby('mech')[['rejects_shared_pi','rejects_shared_k']].mean().round(3))"),
    code("if have('e10_nonstationary_induction.csv'):\n"
         "    display(pd.read_csv(T+'e10_nonstationary_induction.csv').groupby('shift')[\n"
         "        ['true_V_pi_base','true_V_pi_shifted','V_pi']].median().round(4))"),
    md("## 3. The sampling-design result (U5)"),
    code("pd.read_csv(T+'e11_optimal_depth.csv').round(4)"),
    code("if have('e11_equal_budget.csv'):\n"
         "    e = pd.read_csv(T+'e11_equal_budget.csv')\n"
         "    cols = [c for c in ['true_V_pi','V_pi','true_V_logk','V_logk','n_clones_k_identified'] if c in e.columns]\n"
         "    display(e.groupby(['design','C','cells_per_clone','depth'])[cols].median().round(4))"),
    md("### Identifiability windows of two real recording technologies"),
    code("s = pd.read_csv(T+'m3_depth_spectra_summary.csv')\n"
         "s.groupby(['system','dataset'])[['n_leaves','depth_max','k_lo_identifiable',\n"
         "                                 'k_hi_identifiable','window_decades']].median().round(3)"),
    md("## 4. Simeonov re-analysis with the per-clone estimator"),
    code("if have('m3_real_variance_components.csv'):\n"
         "    display(pd.read_csv(T+'m3_real_variance_components.csv').round(4))"),
    code("if have('m3_real_perclone.csv'):\n"
         "    p = pd.read_csv(T+'m3_real_perclone.csv')\n"
         "    display(p[['mouse','clone','n','pi','pi_lo','pi_hi','k','k_at_bound']].round(3).head(30))\n"
         "    print('clones with an identified rate:', int((~p.k_at_bound).sum()), 'of', len(p))"),
    md("### Held-out metastatic site: prediction ladder"),
    code("if have('m3_real_heldout_ladder.csv'):\n"
         "    h = pd.read_csv(T+'m3_real_heldout_ladder.csv')\n"
         "    ec = [c for c in h.columns if c.startswith('err_')]\n"
         "    display(h.groupby('mouse')[ec].mean().round(4))\n"
         "    for m in sorted(h.mouse.unique()):\n"
         "        g = h[(h.mouse==m) & (h.n < h[h.mouse==m].n.max())]\n"
         "        if len(g): print(f'mouse {m} excluding dominant clone:', g[ec].mean().round(4).to_dict())"),
    md("### Leave-one-clone-out stability of the occupancy component"),
    code("if have('m3_real_leave_one_clone_out.csv'):\n"
         "    l = pd.read_csv(T+'m3_real_leave_one_clone_out.csv')\n"
         "    display(l.groupby('mouse').V_pi.describe()[['count','min','50%','max']].round(4))"),
    md("### PATH-style baseline: one autocorrelation number per site"),
    code("if have('m3_real_path_baseline.csv'):\n"
         "    display(pd.read_csv(T+'m3_real_path_baseline.csv').round(4))"),
    md("### Independent-gene validation (held-out programmes, adjusted for clone size, depth, mito, compartment)"),
    code("if have('m3_real_gene_validation.csv'):\n"
         "    g = pd.read_csv(T+'m3_real_gene_validation.csv')\n"
         "    display(g[['mouse','prog','n_clones','rho_pi','p_pi','rho_pi_adj','p_pi_adj']].round(3))"),
]
nb = nbf.v4.new_notebook(cells=cells)
NotebookClient(nb, timeout=1800, kernel_name="python3",
               resources={"metadata": {"path": os.path.join(ROOT, "notebooks")}}).execute()
nbf.write(nb, os.path.join(ROOT, "notebooks/milestone3_results.ipynb"))
print("notebook ok")
