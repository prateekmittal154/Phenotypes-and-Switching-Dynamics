"""Build train/test splits for the held-out prediction benchmark (estimand group D).

Within each clone, 30% of cells are withheld. The training directory contains the tree
pruned to the training cells plus their states, so the existing adapters can be run on it
unchanged. The withheld cells and their true values are stored alongside for scoring.

This is the only comparison in which all four methods appear on one axis, and it is
legitimate because the quantity scored is a prediction about data, not a parameter.

Usage: python -I make_holdout.py <src_sim_dir> <out_dir> <scenarios_csv> <frac> <seed>
"""
import os, sys, csv, glob, json, shutil
import numpy as np

SRC, OUT = sys.argv[1], sys.argv[2]
SCENS = sys.argv[3].split(",")
FRAC = float(sys.argv[4])
SEED = int(sys.argv[5])


def read_newick(path):
    s = open(path).read().strip().rstrip(";")
    parent, elen, labels = [], [], {}
    stack, i = [], 0

    def new_node(p):
        parent.append(p); elen.append(0.0); return len(parent) - 1

    def plen(j):
        if j < len(s) and s[j] == ":":
            k = j + 1
            while k < len(s) and (s[k].isdigit() or s[k] in ".eE+-"):
                k += 1
            return float(s[j + 1:k]), k
        return 0.0, j

    while i < len(s):
        c = s[i]
        if c == "(":
            v = new_node(stack[-1] if stack else -1); stack.append(v); i += 1
        elif c == ",":
            i += 1
        elif c == ")":
            v = stack.pop(); i += 1
            j = i
            while j < len(s) and s[j] not in ":,)":
                j += 1
            i = j
            L, i = plen(i); elen[v] = L
        else:
            j = i
            while j < len(s) and s[j] not in ":,()":
                j += 1
            name = s[i:j]; i = j
            L, i = plen(i)
            v = new_node(stack[-1] if stack else -1)
            elen[v] = L; labels[name] = v
    return np.array(parent, int), np.array(elen, float), labels


def prune(parent, elen, keep_nodes):
    needed = set(keep_nodes)
    for v in keep_nodes:
        x = v
        while parent[x] >= 0:
            x = parent[x]; needed.add(x)
    children = {}
    for v in sorted(needed):
        p = parent[v]
        if p >= 0 and p in needed:
            children.setdefault(p, []).append(v)
    ks = set(keep_nodes)
    retained = {v for v in needed if v in ks or len(children.get(v, [])) >= 2}
    root = min(needed); retained.add(root)
    order = sorted(retained); idx = {v: i for i, v in enumerate(order)}
    npar = np.full(len(order), -1, int); nlen = np.zeros(len(order))
    for v in order:
        if v == root:
            continue
        x, L = parent[v], elen[v]
        while x not in retained:
            L += elen[x]; x = parent[x]
        npar[idx[v]] = idx[x]; nlen[idx[v]] = L
    return npar, nlen, {v: idx[v] for v in keep_nodes if v in retained}


def to_newick(parent, elen, labels):
    kids = {}
    for i in range(len(parent)):
        if parent[i] >= 0:
            kids.setdefault(parent[i], []).append(i)
    root = int(np.where(np.asarray(parent) < 0)[0][0])

    def rec(v):
        if v not in kids:
            return f"{labels[v]}:{max(elen[v], 1e-6):.6f}"
        bl = 0.0 if v == root else max(elen[v], 1e-6)
        return "(" + ",".join(rec(c) for c in kids[v]) + f"):{bl:.6f}"
    return rec(root) + ";"


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)
    n_made = 0
    for scen in SCENS:
        for d in sorted(glob.glob(os.path.join(SRC, f"{scen}_rep*"))):
            rows = list(csv.DictReader(open(os.path.join(d, "states.csv"))))
            parent, elen, labels = read_newick(os.path.join(d, "tree_true.nwk"))
            train, test = [], []
            for cl in sorted({r["clone"] for r in rows}):
                cells = [r for r in rows if r["clone"] == cl]
                idx = rng.permutation(len(cells))
                ncut = max(1, int(round(FRAC * len(cells))))
                test += [cells[i] for i in idx[:ncut]]
                train += [cells[i] for i in idx[ncut:]]
            if len(train) < 24:
                continue
            keep = [labels[r["cell"]] for r in train if r["cell"] in labels]
            npar, nlen, pos = prune(parent, elen, keep)
            lab = [None] * len(npar)
            for r in train:
                v = labels.get(r["cell"])
                if v is not None and v in pos:
                    lab[pos[v]] = r["cell"]

            od = os.path.join(OUT, os.path.basename(d))
            os.makedirs(od, exist_ok=True)
            open(os.path.join(od, "tree_true.nwk"), "w").write(to_newick(npar, nlen, lab) + "\n")
            open(os.path.join(od, "tree_unit.nwk"), "w").write(
                to_newick(npar, np.ones(len(npar)), lab) + "\n")
            with open(os.path.join(od, "states.csv"), "w", newline="") as f:
                w = csv.writer(f); w.writerow(["cell", "clone", "state", "trait"])
                for r in train:
                    w.writerow([r["cell"], r["clone"], r["state"], r["trait"]])
            with open(os.path.join(od, "heldout.csv"), "w", newline="") as f:
                w = csv.writer(f); w.writerow(["cell", "clone", "state", "trait"])
                for r in test:
                    w.writerow([r["cell"], r["clone"], r["state"], r["trait"]])
            shutil.copy(os.path.join(d, "truth.json"), os.path.join(od, "truth.json"))
            n_made += 1
    print(f"wrote {n_made} held-out replicates to {OUT}")


if __name__ == "__main__":
    main()
