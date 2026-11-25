"""Build and execute notebooks/milestone4_results.ipynb."""
import os, nbformat as nbf
from nbclient import NotebookClient
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = [
    md("# Study D, Milestone 4\n"
       "Dataset expansion, audited Fisher information, full-tree dependence, calibration, the capacity law "
       "verified against real tree geometry, and the Quinn replication. See `reports/milestone4_decision.md`."),
    code("import json, os, sys, numpy as np, pandas as pd\nsys.path.insert(0, '..')\n"
         "pd.set_option('display.width', 220)\nT = '../results/tables/'\n"
         "have = lambda f: os.path.exists(T + f)"),
    md("## 1. Fisher-information audit: the Milestone 3 optimum restated"),
    code("o = pd.read_csv(T+'e12_optimum_audit.csv')\n"
         "display(o.round(4))\n"
         "print('kL* (pi known):  mean %.3f  range %.3f-%.3f' % (o.kL_opt_pi_known.mean(),\n"
         "      o.kL_opt_pi_known.min(), o.kL_opt_pi_known.max()))\n"
         "print('kL* (pi joint):  mean %.3f  range %.3f-%.3f' % (o.kL_opt_pi_joint.mean(),\n"
         "      o.kL_opt_pi_joint.min(), o.kL_opt_pi_joint.max()))\n"
         "print('worst-case information cost of estimating pi: %.1f%%' %\n"
         "      (100*(1-o.ratio_joint_over_known.min())))"),
    code("w = pd.read_csv(T+'e12_optimum_width.csv')\n"
         "print('kL band retaining >=50%% of peak: %.2f - %.2f (%.1f-fold)' %\n"
         "      (w.kL_lo_50pct.median(), w.kL_hi_50pct.median(), w.fold_width.median()))"),
    md("## 2. Full-tree information versus the naive independent-pair sum"),
    code("if have('e13_tree_vs_pair_mc.csv'):\n"
         "    d = pd.read_csv(T+'e13_tree_vs_pair_mc.csv')\n"
         "    display(d.groupby('depth').inflation.describe()[['count','min','50%','max']].round(2))\n"
         "    display(d.groupby('k').inflation.median().round(2).rename('median inflation'))\n"
         "    v = d.dropna(subset=['I_exact_logk'])\n"
         "    print('MC vs exact enumeration: max relative error %.4f over %d settings' %\n"
         "          (((v.I_tree_logk-v.I_exact_logk).abs()/v.I_exact_logk).max(), len(v)))"),
    md("## 3. Calibration (300 replicates per scenario, Wilson intervals)"),
    code("if have('e14_calibration_summary.csv'):\n"
         "    s = pd.read_csv(T+'e14_calibration_summary.csv')\n"
         "    print('replicates per scenario:', int(s.n.iloc[0]))\n"
         "    display(s.pivot(index='scen', columns='test', values='rate').round(3))\n"
         "    lo = s.pivot(index='scen', columns='test', values='lo')\n"
         "    hi = s.pivot(index='scen', columns='test', values='hi')\n"
         "    for sc in [x for x in ['H0_pi','H0_k','H0_k_frozen','MIS'] if x in lo.index]:\n"
         "        print(sc, {t: '%.3f [%.3f, %.3f]' % (s.pivot(index='scen',columns='test',values='rate').loc[sc,t],\n"
         "              lo.loc[sc,t], hi.loc[sc,t]) for t in ['chi2_pi','chi2_k','cover_V_pi'] if t in lo})"),
    code("if have('e14_thresholds.csv'):\n    display(pd.read_csv(T+'e14_thresholds.csv').round(2))"),
    md("## 4. The capacity law, verified against real tree geometry"),
    code("if have('e16_capacity_vs_prediction.csv'):\n"
         "    e = pd.read_csv(T+'e16_capacity_vs_prediction.csv')\n"
         "    display(e.round(3))\n"
         "    print('median |empirical - predicted| decades: %.3f' % e.abs_error.median())"),
    code("if have('e15b_capacity_designs.csv'):\n"
         "    d = pd.read_csv(T+'e15b_capacity_designs.csv')\n"
         "    display(d.groupby(['S','D']).decades_covered.first().unstack().round(2))\n"
         "    print('Pareto-optimal designs:', int(d.pareto.sum()), 'of', len(d))"),
    md("## 5. Quinn independent system: four animals"),
    code("for f, lab in [('m4_quinn_summary.csv','cohort'),\n"
         "               ('m4_quinn_variance_components.csv','all sites, cell-cycle state'),\n"
         "               ('m4_quinn_variance_components_pc1.csv','all sites, PC1 sensitivity'),\n"
         "               ('m4_quinn_primary_only.csv','primary site only, depth-adjusted')]:\n"
         "    if have(f):\n"
         "        print('==', lab); display(pd.read_csv(T+f).round(4))"),
    md("### Negative control: within-clone shifts against quality covariates"),
    code("if have('m4_quinn_technical_control.csv'):\n"
         "    t = pd.read_csv(T+'m4_quinn_technical_control.csv')\n"
         "    display(t.pivot_table(index=['mouse','site'], columns='feature',\n"
         "                          values='within_clone_shift_sd').round(3))"),
]
nb = nbf.v4.new_notebook(cells=cells)
NotebookClient(nb, timeout=1800, kernel_name="python3",
               resources={"metadata": {"path": os.path.join(ROOT, "notebooks")}}).execute()
nbf.write(nb, os.path.join(ROOT, "notebooks/milestone4_results.ipynb"))
print("notebook ok")
