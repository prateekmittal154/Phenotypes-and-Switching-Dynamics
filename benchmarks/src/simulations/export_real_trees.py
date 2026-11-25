"""Export the real Simeonov and Quinn lineage trees as Newick plus state tables.

Mirrors the loaders in the canonical Study D repository
(experiments/reconstruction/run_e22_realboot.py) so the competing methods see exactly the
trees and state definitions Study D used. Both the binary state and the underlying
continuous score are written, because SCOUT needs a continuous trait.

Branch lengths: the published trees carry no usable branch times, so edges are written
with length 1. That is the same convention SCOUT applies internally
(R/dropout_validation.R:470-472) and the one PASTRI implies by using topological depth.

Usage: python -I export_real_trees.py <studyd_src> <outdir>
"""
import os, sys, re, glob, csv
import numpy as np
import pandas as pd

S = sys.argv[1]
OUT = sys.argv[2]
sys.path.insert(0, S)

MIN_CELLS, MAX_CELLS = 25, 400

PROC = os.path.join(S, "data", "simeonov")
CTAB = os.path.join(S, "results", "tables")
QTREES = os.path.join(S, "data", "quinn_trees")
QSTATES = os.path.join(S, "data", "quinn_states")
GSM = {"5k": "GSM4905334_trees_5k", "30k": "GSM4905336_trees_30k",
       "10k": "GSM4905338_trees_10k", "100k": "GSM4905340_trees_100k"}


def newick_from_parent(parent, labels):
    n = len(parent)
    kids = {}
    for i in range(n):
        if parent[i] >= 0:
            kids.setdefault(parent[i], []).append(i)
    roots = [i for i in range(n) if parent[i] < 0]
    if len(roots) != 1:
        raise ValueError(f"newick_from_parent: {len(roots)} roots, expected 1")
    root = roots[0]

    def rec(v):
        if v not in kids:
            lab = labels[v] if labels[v] else f"n{v}"
            return f"{lab}:1.0"
        return "(" + ",".join(rec(c) for c in kids[v]) + "):1.0"
    return rec(root) + ";"


def prune_to_labelled(parent, labels):
    """Keep only labelled leaves and the internal structure joining them.

    The root is the node whose parent is -1. It must NOT be taken as min(needed): the
    Simeonov nodes are indexed by sorted name, so the lowest index is the alphabetically
    first node, not the root. Getting this wrong collapsed whole clones to a single tip.
    """
    keep_leaves = [i for i, l in enumerate(labels) if l]
    needed = set(keep_leaves)
    for v in keep_leaves:
        x = v
        while parent[x] >= 0:
            x = parent[x]
            needed.add(x)
    children = {}
    for v in sorted(needed):
        p = parent[v]
        if p >= 0 and p in needed:
            children.setdefault(p, []).append(v)
    keepset = set(keep_leaves)
    retained = {v for v in needed if v in keepset or len(children.get(v, [])) >= 2}
    roots = [v for v in needed if parent[v] < 0]
    if len(roots) != 1:
        raise ValueError(f"expected exactly one root, found {len(roots)}")
    root = roots[0]
    retained.add(root)
    order = sorted(retained)
    idx = {v: i for i, v in enumerate(order)}
    npar = np.full(len(order), -1, int)
    for v in order:
        if v == root:
            continue
        x = parent[v]
        while x not in retained:
            x = parent[x]
        npar[idx[v]] = idx[x]
    nlab = [labels[v] if v in keepset else "" for v in order]
    return npar, nlab


def write_one(tag, parent, labels, cells, binary, cont):
    if len(cells) < 2 or len(parent) < 3:
        raise ValueError(f"{tag}: degenerate tree with {len(cells)} labelled tips")
    d = os.path.join(OUT, tag)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "tree_unit.nwk"), "w") as f:
        f.write(newick_from_parent(parent, labels) + "\n")
    with open(os.path.join(d, "states.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cell", "clone", "state", "trait"])
        for c in cells:
            w.writerow([c, tag.split("__")[1], int(binary[c]), f"{float(cont[c]):.6f}"])
    return len(cells)


def do_simeonov():
    rows = []
    for m in (1, 2):
        d = pd.read_csv(os.path.join(CTAB, f"m2_real_M{m}_cells.csv.gz"))
        el = pd.read_csv(os.path.join(PROC, f"M{m}/output-files/pdac_mouse{m}_full_edgelist.txt"),
                         sep="\t")
        el["V1"] = el.V1.str.strip(); el["V2"] = el.V2.str.strip()
        binary = dict(zip(d.leaf.astype(str), d.state2))
        cont = dict(zip(d.leaf.astype(str), d["e"]))
        for c, g in el.groupby("cloneID"):
            par = dict(zip(g.V2, g.V1))
            nodes = sorted(set(g.V1) | set(g.V2))
            obs = [n for n in nodes if n in binary]
            if not (MIN_CELLS <= len(obs) <= MAX_CELLS):
                continue
            keep = set(obs)
            for n in obs:
                x = n
                while x in par:
                    x = par[x]; keep.add(x)
            keep = sorted(keep); idx = {n: i for i, n in enumerate(keep)}
            parent = np.array([idx[par[n]] if (n in par and par[n] in idx) else -1 for n in keep])
            if int((parent < 0).sum()) != 1:
                continue
            labels = [n if n in binary else "" for n in keep]
            parent, labels = prune_to_labelled(parent, labels)
            cells = [l for l in labels if l]
            tag = f"simeonov_M{m}__{int(c)}"
            n = write_one(tag, parent, labels, cells, binary, cont)
            rows.append(dict(dataset=f"Simeonov M{m}", animal=f"M{m}", clone=int(c),
                             n_cells=n, tag=tag))
    return rows


def do_quinn():
    from experiments.quinn_analysis.run_quinn import PRIMARY, parse_newick
    rows = []
    for mouse, gsm in GSM.items():
        f = os.path.join(QSTATES, f"quinn_cells_{mouse}.csv.gz")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        d["site"] = d.sampleID.astype(str).str.replace("^30k", "", regex=True)
        d = d.assign(state_use=d.state2_cc)
        for tf in sorted(glob.glob(os.path.join(QTREES, gsm, "trees", "*hybrid_priors*.txt"))):
            mm = re.search(r"lg(\d+)", os.path.basename(tf))
            if mm is None:
                continue
            lg = int(mm.group(1)); dl = d[d.LineageGroup == lg]
            if len(dl) < MIN_CELLS:
                continue
            s1b = dict(zip(dl.cellBC.astype(str), dl.state_use))
            s1c = dict(zip(dl.cellBC.astype(str), dl.cc_score))
            s2b = dict(zip(dl.site.astype(str) + "." + dl.cellBC.astype(str), dl.state_use))
            s2c = dict(zip(dl.site.astype(str) + "." + dl.cellBC.astype(str), dl.cc_score))
            par, nm = parse_newick(open(tf).read())
            best = None
            for sb, sc in ((s1b, s1c), (s2b, s2c)):
                hit = sum(1 for i, name in nm.items() if name in sb)
                if best is None or hit > best[0]:
                    best = (hit, sb, sc)
            hit, binary, cont = best
            if not (MIN_CELLS <= hit <= MAX_CELLS):
                continue
            labels = ["" for _ in range(len(par))]
            for i, name in nm.items():
                if name in binary:
                    labels[i] = name
            parent, labels = prune_to_labelled(np.asarray(par), labels)
            cells = [l for l in labels if l]
            tag = f"quinn_{mouse}__{lg}"
            n = write_one(tag, parent, labels, cells, binary, cont)
            rows.append(dict(dataset=f"Quinn {mouse}", animal=mouse, clone=lg,
                             n_cells=n, tag=tag))
    return rows


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    rows = do_simeonov() + do_quinn()
    man = pd.DataFrame(rows)
    man.to_csv(os.path.join(OUT, "manifest.csv"), index=False)
    print(man.groupby("dataset").agg(clones=("clone", "size"),
                                     cells=("n_cells", "sum")).to_string())
    print(f"\nwrote {len(rows)} clone trees to {OUT}")
