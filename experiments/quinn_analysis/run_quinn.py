"""Quinn (GSE161363) independent-system analysis: occupancy, switching-rate identifiability, composition.

Four independent mice, each with its own barcode library complexity (5k, 10k, 30k, 100k) and several anatomical
sites. States were derived within this system by `quinn_states.py` (unsupervised PC1 of primary-site cells),
never from the Simeonov EMT definition.

Usage: python -I run_quinn.py <quinn_dir> <states_dir> <tables_dir>
"""
import os, re, sys, glob
import numpy as np, pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.lineage_inference.ctmc import FastTree
from src.lineage_inference.perclone import fit_clone, variance_components, lrt_shared

TREES = {"5k": "GSM4905334_trees_5k", "30k": "GSM4905336_trees_30k",
         "10k": "GSM4905338_trees_10k", "100k": "GSM4905340_trees_100k"}
PRIMARY = ("LL", "RL", "RE", "RW", "NE", "NW", "SE", "SW")
rng = np.random.default_rng(20261202)


def parse_newick(s):
    """Return (parent array, leaf node indices, leaf names). Internal labels and branch lengths are ignored."""
    par = [-1]; cur = 0; stack = [0]
    names = {}
    i = 0
    while i < len(s):
        ch = s[i]
        if ch == "(":
            par.append(stack[-1]); stack.append(len(par) - 1); i += 1
        elif ch == ")":
            stack.pop()
            j = i + 1
            while j < len(s) and s[j] not in ",();":
                j += 1
            i = j
        elif ch in ",;":
            if ch == ",":
                # the next token is a sibling of the current node: attach to the enclosing parent
                pass
            i += 1
        else:
            j = i
            while j < len(s) and s[j] not in ",();":
                j += 1
            nm = s[i:j].strip()
            if nm:
                par.append(stack[-1]); names[len(par) - 1] = nm.split(":")[0]
            i = j
    return np.array(par), names


MAX_LEAVES = 1500     # trees larger than this are skipped: see reports/milestone4_dataset_expansion.md
MAX_TREES = 30        # per mouse, largest first, to bound the likelihood-ratio optimisation


def build(nwk, state_of):
    par, names = parse_newick(nwk)
    ls = np.full(len(par), -1)
    n_match = 0
    for idx, nm in names.items():
        if nm in state_of:
            ls[idx] = state_of[nm]; n_match += 1
    if n_match < 25 or n_match > MAX_LEAVES:
        return None, n_match
    return FastTree(par, ls), n_match


def main(qdir, sdir, tab, state_col="state2_cc"):
    """state_col: 'state2_cc' (cell cycle, the primary definition) or 'state2' (unsupervised PC1, secondary).

    PC1 was found to correlate with mitochondrial fraction (-0.67 to +0.60 across mice) even after excluding
    mitochondrial and ribosomal genes, so it is substantially a cell-quality axis and is reported only as a
    sensitivity analysis. The cell-cycle score uses a fixed external gene set, is biologically interpretable,
    and is not fitted to these data.
    """
    os.makedirs(tab, exist_ok=True)
    rows_pc, rows_vc, rows_comp, rows_sum = [], [], [], []
    for mouse, tdir in TREES.items():
        f = os.path.join(sdir, f"quinn_cells_{mouse}.csv.gz")
        if not os.path.exists(f):
            print("missing states for", mouse); continue
        d = pd.read_csv(f)
        d["site"] = d.sampleID.astype(str).str.replace("^30k", "", regex=True)
        d["is_primary"] = d.site.isin(PRIMARY)
        # the tree leaf names may or may not carry the "<site>." prefix; build both lookups
        if state_col not in d.columns:
            print("missing", state_col, "for", mouse); continue
        trees, names, nmatch = [], [], []
        # three reconstructions exist per lineage group (hybrid/nj/nj_priors); use the authors' primary
        # hybrid-with-priors tree only, otherwise every cell is counted three times
        for tf in sorted(glob.glob(os.path.join(qdir, tdir, "trees", "*hybrid_priors*.txt"))):
            nwk = open(tf).read()
            lg = re.search(r"lg(\d+)", os.path.basename(tf))
            if lg is None:
                continue
            lgi = int(lg.group(1))
            dl = d[d.LineageGroup == lgi]
            if len(dl) < 25:
                continue
            s1 = dict(zip(dl.cellBC.astype(str), dl[state_col]))
            s2 = dict(zip(dl.site.astype(str) + "." + dl.cellBC.astype(str), dl[state_col]))
            t1, n1 = build(nwk, s1); t2, n2 = build(nwk, s2)
            t, n = (t1, n1) if n1 >= n2 else (t2, n2)
            if t is not None:
                trees.append(t); names.append(lgi); nmatch.append(n)
        rows_sum.append(dict(mouse=mouse, state=state_col, n_cells=len(d), n_sites=int(d.site.nunique()),
                             n_lineage_groups=int(d.LineageGroup.nunique()),
                             n_trees_usable=len(trees), cells_in_trees=int(sum(nmatch)),
                             primary_cells=int(d.is_primary.sum())))
        print(mouse, "usable trees", len(trees), "cells in trees", sum(nmatch), flush=True)
        if len(trees) > MAX_TREES:
            keep = np.argsort(-np.array(nmatch))[:MAX_TREES]
            trees = [trees[i] for i in keep]; names = [names[i] for i in keep]
            nmatch = [nmatch[i] for i in keep]
        if len(trees) >= 3:
            fits = [fit_clone(t) for t in trees]
            for c, fi, nm in zip(names, fits, nmatch):
                rows_pc.append(dict(mouse=mouse, state=state_col, clone=c, n_matched=nm, **fi))
            vc = variance_components(fits, min_n=25)
            lp = lrt_shared(trees, fits, "pi", maxiter_mult=25); lk = lrt_shared(trees, fits, "k", maxiter_mult=25)
            rows_vc.append(dict(mouse=mouse, state=state_col, scope="all_sites", lrt_pi_stat=lp["stat"], lrt_pi_p=lp["p_asym"],
                                lrt_k_stat=lk["stat"], lrt_k_p=lk["p_asym"], **vc))
        # composition versus within-clone decomposition, per disseminated site
        base = d[d.is_primary]
        for s in sorted(set(d.site) - set(PRIMARY)):
            sub = d[d.site == s]
            if len(sub) < 40 or len(base) < 40:
                continue
            pc = base.LineageGroup.value_counts(normalize=True)
            qc = sub.LineageGroup.value_counts(normalize=True)
            mu = base.groupby("LineageGroup")[state_col].mean()
            mus = sub.groupby("LineageGroup")[state_col].mean()
            cl = [c for c in set(pc.index) | set(qc.index) if c in mu.index]
            comp = sum(qc.get(c, 0) * mu[c] for c in cl) - sum(pc.get(c, 0) * mu[c] for c in cl)
            within = sum(qc.get(c, 0) * (mus.get(c, np.nan) - mu[c]) for c in cl
                         if c in mus.index and np.isfinite(mus.get(c, np.nan)))
            rows_comp.append(dict(mouse=mouse, state=state_col, site=s, n_site=len(sub), n_primary=len(base),
                                  total=float(sub[state_col].mean() - base[state_col].mean()),
                                  composition=float(comp), within_clone=float(within),
                                  n_shared_clones=len(cl)))
    tag = "" if state_col == "state2_cc" else "_pc1"
    pd.DataFrame(rows_sum).to_csv(os.path.join(tab, f"m4_quinn_summary{tag}.csv"), index=False)
    pd.DataFrame(rows_pc).to_csv(os.path.join(tab, f"m4_quinn_perclone{tag}.csv"), index=False)
    pd.DataFrame(rows_vc).to_csv(os.path.join(tab, f"m4_quinn_variance_components{tag}.csv"), index=False)
    pd.DataFrame(rows_comp).to_csv(os.path.join(tab, f"m4_quinn_composition{tag}.csv"), index=False)
    for nm, fr in [("summary", rows_sum), ("variance", rows_vc), ("composition", rows_comp)]:
        if fr:
            print("==", nm); print(pd.DataFrame(fr).round(4).to_string())
    print("done")


if __name__ == "__main__":
    main(*sys.argv[1:4], state_col=(sys.argv[4] if len(sys.argv) > 4 else "state2_cc"))
