"""Cross-technology comparison of MRCA-depth spectra and the identifiability window (Proposition U4).

For each lineage-recording system, the distribution of shared-ancestor depth over within-clone leaf pairs sets
the window 1/(2 L_max) <~ k <~ 1/(2 L_min) inside which a switching rate is identifiable. This script computes
that spectrum for macsGESTALT (Simeonov) and Cassiopeia (Quinn) and reports the implied window in units of
branch length.

Usage: python -I depth_spectra.py <tables_dir> <proc_dir> <quinn_dir>
"""
import os, re, sys, glob
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
rng = np.random.default_rng(20261131)


def newick_leaf_depths(s):
    """Depth of each leaf (number of enclosing parentheses) and pairwise MRCA depths, from a newick string."""
    depth, leaves, buf, cur = 0, [], "", 0
    i = 0
    paths = []
    stack = []
    while i < len(s):
        ch = s[i]
        if ch == "(":
            cur += 1; stack.append(0); i += 1
        elif ch == ")":
            cur -= 1
            if stack:
                stack.pop()
            # skip the internal node label up to the next , ) or ;
            j = i + 1
            while j < len(s) and s[j] not in ",();":
                j += 1
            i = j
        elif ch == ",":
            i += 1
        else:
            j = i
            while j < len(s) and s[j] not in ",();":
                j += 1
            name = s[i:j].strip()
            if name:
                leaves.append(name); paths.append(cur)
            i = j
    return leaves, np.array(paths)


def pair_mrca_depths_from_nested(s, max_leaves=400):
    """Sample leaf pairs and return their MRCA depth, by tracking the parenthesis path of each leaf."""
    path, cur, out_names, out_paths = [], 0, [], []
    counter = 0
    i = 0
    while i < len(s):
        ch = s[i]
        if ch == "(":
            counter += 1; path.append(counter); i += 1
        elif ch == ")":
            if path:
                path.pop()
            j = i + 1
            while j < len(s) and s[j] not in ",();":
                j += 1
            i = j
        elif ch in ",;":
            i += 1
        else:
            j = i
            while j < len(s) and s[j] not in ",();":
                j += 1
            nm = s[i:j].strip()
            if nm:
                out_names.append(nm); out_paths.append(tuple(path))
            i = j
    if len(out_names) > max_leaves:
        k = rng.choice(len(out_names), max_leaves, replace=False)
        out_names = [out_names[t] for t in k]; out_paths = [out_paths[t] for t in k]
    d = []
    for a in range(len(out_paths)):
        for b in range(a + 1, len(out_paths)):
            x, y = out_paths[a], out_paths[b]
            sh = 0
            for u, v in zip(x, y):
                if u != v:
                    break
                sh += 1
            d.append(sh)
    return np.array(d), len(out_names)


def main(tab, proc, quinn):
    rows, spec = [], []
    # --- macsGESTALT (Simeonov): depths come from the published edgelist
    for m in (1, 2):
        el = pd.read_csv(os.path.join(proc, f"M{m}/output-files/pdac_mouse{m}_full_edgelist.txt"), sep="\t")
        el["V1"] = el.V1.str.strip(); el["V2"] = el.V2.str.strip()
        par = dict(zip(el.V2, el.V1))
        leaves = sorted(set(el.V2) - set(el.V1))
        paths = {}
        for l in leaves:
            p, n = [], l
            while n in par:
                n = par[n]
                if n != "Cancer_root":
                    p.append(n)
            paths[l] = tuple(reversed(p))
        cl = dict(zip(el.V2, el.cloneID))
        for c, g in pd.Series(cl).groupby(lambda x: cl[x]):
            ls = [l for l in g.index if l in paths]
            if len(ls) < 25:
                continue
            if len(ls) > 400:
                ls = [ls[t] for t in rng.choice(len(ls), 400, replace=False)]
            d = []
            for a in range(len(ls)):
                for b in range(a + 1, len(ls)):
                    x, y = paths[ls[a]], paths[ls[b]]
                    sh = 0
                    for u, v in zip(x, y):
                        if u != v:
                            break
                        sh += 1
                    d.append(sh)
            d = np.array(d)
            if len(d) < 50:
                continue
            rows.append(dict(system="macsGESTALT", dataset=f"mouse{m}", clone=str(c), n_leaves=len(ls),
                             depth_min=int(d.min()), depth_max=int(d.max()),
                             depth_p10=float(np.percentile(d, 10)), depth_p90=float(np.percentile(d, 90)),
                             depth_distinct=int(len(np.unique(d))), n_pairs=len(d)))
            for h in np.unique(d):
                spec.append(dict(system="macsGESTALT", dataset=f"mouse{m}", clone=str(c), depth=int(h),
                                 frac=float((d == h).mean())))
    # --- Cassiopeia (Quinn): newick trees
    for f in sorted(glob.glob(os.path.join(quinn, "*trees*", "trees", "*.txt"))):
        s = open(f).read()
        d, nl = pair_mrca_depths_from_nested(s)
        if len(d) < 50:
            continue
        tag = re.sub(r"_tree.*", "", os.path.basename(f))
        rows.append(dict(system="Cassiopeia", dataset=tag.split("_")[0], clone=tag, n_leaves=nl,
                         depth_min=int(d.min()), depth_max=int(d.max()),
                         depth_p10=float(np.percentile(d, 10)), depth_p90=float(np.percentile(d, 90)),
                         depth_distinct=int(len(np.unique(d))), n_pairs=len(d)))
        for h in np.unique(d):
            spec.append(dict(system="Cassiopeia", dataset=tag.split("_")[0], clone=tag, depth=int(h),
                             frac=float((d == h).mean())))
    r = pd.DataFrame(rows)
    if len(r):
        # identifiability window in branch-length units, using L = depth_max - depth as lineage distance
        r["L_max"] = r.depth_max - r.depth_p10
        r["L_min"] = np.maximum(r.depth_max - r.depth_p90, 0.5)
        r["k_lo_identifiable"] = 1 / (2 * r.L_max.clip(lower=0.5))
        r["k_hi_identifiable"] = 1 / (2 * r.L_min)
        r["window_decades"] = np.log10(r.k_hi_identifiable / r.k_lo_identifiable)
    r.to_csv(os.path.join(tab, "m3_depth_spectra_summary.csv"), index=False)
    pd.DataFrame(spec).to_csv(os.path.join(tab, "m3_depth_spectra.csv"), index=False)
    if len(r):
        print(r.groupby(["system", "dataset"])[
            ["n_leaves", "depth_max", "depth_distinct", "k_lo_identifiable", "k_hi_identifiable",
             "window_decades"]].median().round(3).to_string())
    print("done")


if __name__ == "__main__":
    main(*sys.argv[1:4])
