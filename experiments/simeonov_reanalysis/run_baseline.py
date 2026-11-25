"""PATH-style phylogenetic autocorrelation as the established-method baseline, on the real Simeonov data.

A single depth-weighted autocorrelation per mouse-site. It is reported alongside the per-clone decomposition to
show what a one-number summary can and cannot distinguish.
"""
import os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.lineage_inference.perclone import moran_path

rng = np.random.default_rng(20261132)


def main(tab, proc):
    rows = []
    for m in (1, 2):
        d = pd.read_csv(os.path.join(tab, f"m2_real_M{m}_cells.csv.gz"))
        el = pd.read_csv(os.path.join(proc, f"M{m}/output-files/pdac_mouse{m}_full_edgelist.txt"), sep="\t")
        el["V1"] = el.V1.str.strip(); el["V2"] = el.V2.str.strip()
        par = dict(zip(el.V2, el.V1))

        def path(leaf):
            p, n = [], leaf
            while n in par:
                n = par[n]
                if n != "Cancer_root":
                    p.append(n)
            return tuple(reversed(p))

        for s in ["PT", "Met", "Liver"]:
            sub = d[d.site == s]
            if len(sub) < 150:
                continue
            k = rng.choice(len(sub), min(len(sub), 900), replace=False)
            sm = sub.iloc[k]
            paths = [path(l) for l in sm.leaf]
            val = moran_path(paths, sm.state2.values, sm.cloneID.values)
            rows.append(dict(mouse=m, site=s, n_cells=len(sm), path_autocorrelation=val,
                             n_clones=int(sm.cloneID.nunique())))
            print(m, s, round(val, 4), flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(tab, "m3_real_path_baseline.csv"), index=False)
    print("done")


if __name__ == "__main__":
    main(*sys.argv[1:3])
