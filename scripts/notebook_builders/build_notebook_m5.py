"""Build and execute notebooks/milestone5_results.ipynb."""
import os, nbformat as nbf
from nbclient import NotebookClient
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = [
    md("# Study D, Milestone 5\n"
       "Separating biological sampling loss, recorder loss and tree-reconstruction error. See "
       "`reports/milestone5_decision.md`."),
    code("import os, sys, numpy as np, pandas as pd\nsys.path.insert(0, '..')\n"
         "pd.set_option('display.width', 230)\nT = '../results/tables/'\n"
         "have = lambda f: os.path.exists(T + f)"),
    md("## 1. Where switching information is lost (E17, 480 replicates)\n"
       "L1 true tree with true branch lengths; L2 same topology, unit edges; L3 reconstructed."),
    code("d = pd.read_csv(T+'e17_hierarchy.csv')\n"
         "print('replicates:', len(d))\n"
         "display(pd.read_csv(T+'e17_paired_ratios.csv').round(3))"),
    md("### Star-likeness is already present in the TRUE sampled genealogy"),
    code("cols = [c for c in d.columns if c.startswith('close_frac')]\n"
         "display(d.groupby(['recorder','genealogy'])[cols].median().round(4))"),
    md("### Topology recovery, and what recorder capacity buys"),
    code("rc = [c for c in d.columns if c.startswith('rf_recall')]\n"
         "display(d.groupby(['recorder','genealogy'])[rc].median().round(3))\n"
         "display(d.groupby('recorder')[['saturation','missing','n_informative_chars']].median().round(3))"),
    md("### Ensemble versus single trees"),
    code("for ref in ['k_L3_greedy','k_L3_nj','k_L3_upgma']:\n"
         "    b = (np.abs(np.log(d.k_L3_ens/d.k_true)) < np.abs(np.log(d[ref]/d.k_true))).mean()\n"
         "    print(f'ensemble closer to truth than {ref}: {b:.3f}')\n"
         "es = np.abs(np.log(d[['k_L3_greedy','k_L3_nj','k_L3_upgma']].values/d.k_true.values[:,None])).mean(1)\n"
         "ee = np.abs(np.log(d.k_L3_ens/d.k_true))\n"
         "print('mean |log error|: single %.3f, ensemble %.3f' % (es.mean(), ee.mean()))"),
    md("## 2. Real data: three published reconstructions of the same cells (Quinn 5k)"),
    code("if have('m5_quinn_variance_by_reconstruction.csv'):\n"
         "    display(pd.read_csv(T+'m5_quinn_variance_by_reconstruction.csv').round(4))\n"
         "if have('m5_quinn_reconstruction_geometry.csv'):\n"
         "    g = pd.read_csv(T+'m5_quinn_reconstruction_geometry.csv')\n"
         "    display(g.groupby(['mouse','variant'])[['height','mrca_norm','close_frac','L_min','L_max']]\n"
         "             .median().round(3))"),
    md("## 3. Real data: Simeonov, reconstructed from the authors' allele records"),
    code("if have('m5_simeonov_variance_by_reconstruction.csv'):\n"
         "    display(pd.read_csv(T+'m5_simeonov_variance_by_reconstruction.csv').round(4))\n"
         "if have('m5_simeonov_alt_reconstruction_fits.csv'):\n"
         "    r = pd.read_csv(T+'m5_simeonov_alt_reconstruction_fits.csv')\n"
         "    display(r.groupby(['mouse','method']).agg(clones=('clone','size'), pi=('pi','median'),\n"
         "            k=('k','median'), at_bound=('k_at_bound','mean')).round(3))"),
    md("## 4. The Milestone 4 calibration, completed (300 replicates)"),
    code("if have('e14_calibration_summary.csv'):\n"
         "    s = pd.read_csv(T+'e14_calibration_summary.csv')\n"
         "    print('replicates per scenario:', int(s.n.iloc[0]))\n"
         "    display(s.pivot(index='scen', columns='test', values='rate').round(3))\n"
         "    display(pd.read_csv(T+'e14_thresholds.csv').round(2))"),
]
nb = nbf.v4.new_notebook(cells=cells)
NotebookClient(nb, timeout=1800, kernel_name="python3",
               resources={"metadata": {"path": os.path.join(ROOT, "notebooks")}}).execute()
nbf.write(nb, os.path.join(ROOT, "notebooks/milestone5_results.ipynb"))
print("notebook ok")
