"""Run the Study D estimator on the same replicates the competitors see.

Study D fits (pi_c, k_c) per clone by maximum likelihood with Felsenstein pruning, and
reports profile-likelihood intervals. Both the true-branch-length and unit-edge
conventions are supported so the clock comparison is like-for-like.

Usage: python -I run_studyd.py <studyd_src> <list_file> <tree_variant> <out_csv>
"""
import os, sys, json, csv, time
import numpy as np

STUDYD = sys.argv[1]
LIST = sys.argv[2]
TREEV = sys.argv[3]
OUT = sys.argv[4]

sys.path.insert(0, STUDYD)
from src.lineage_inference.ctmc import VarTree, fit_var
from src.lineage_inference.perclone import fit_clone


def read_newick_topology(path):
    """Minimal Newick parser -> (parent array, edge lengths, leaf label -> node index)."""
    s = open(path).read().strip().rstrip(";")
    parent, elen, labels = [], [], {}
    stack = []
    i = 0
    cur_root = None

    def new_node(p):
        parent.append(p)
        elen.append(0.0)
        return len(parent) - 1

    def parse_len(j):
        if j < len(s) and s[j] == ":":
            k = j + 1
            while k < len(s) and (s[k].isdigit() or s[k] in ".eE+-"):
                k += 1
            return float(s[j + 1:k]), k
        return 0.0, j

    while i < len(s):
        c = s[i]
        if c == "(":
            p = stack[-1] if stack else -1
            v = new_node(p)
            if p == -1:
                cur_root = v
            stack.append(v)
            i += 1
        elif c == ",":
            i += 1
        elif c == ")":
            v = stack.pop()
            i += 1
            # optional internal label
            j = i
            while j < len(s) and s[j] not in ":,)":
                j += 1
            i = j
            L, i = parse_len(i)
            elen[v] = L
            if not stack:
                cur_root = v
        else:
            j = i
            while j < len(s) and s[j] not in ":,()":
                j += 1
            name = s[i:j]
            i = j
            L, i = parse_len(i)
            p = stack[-1] if stack else -1
            v = new_node(p)
            elen[v] = L
            labels[name] = v
    _ = cur_root
    return np.array(parent, int), np.array(elen, float), labels


def subtree_for(parent, elen, node_idx):
    """Induced subtree on the given leaf nodes, collapsing unary chains."""
    needed = set(node_idx)
    for v in node_idx:
        x = v
        while parent[x] >= 0:
            x = parent[x]
            needed.add(x)
    children = {}
    for v in sorted(needed):
        p = parent[v]
        if p >= 0 and p in needed:
            children.setdefault(p, []).append(v)
    keepset = set(node_idx)
    retained = {v for v in needed if v in keepset or len(children.get(v, [])) >= 2}
    root = min(needed)
    retained.add(root)
    order = sorted(retained)
    idx = {v: i for i, v in enumerate(order)}
    npar = np.full(len(order), -1, int)
    nlen = np.zeros(len(order), float)
    for v in order:
        if v == root:
            continue
        x = parent[v]
        L = elen[v]
        while x not in retained:
            L += elen[x]
            x = parent[x]
        npar[idx[v]] = idx[x]
        nlen[idx[v]] = L
    leafpos = {v: idx[v] for v in node_idx if v in retained}
    return npar, nlen, leafpos


def main():
    reps = [l.strip() for l in open(LIST) if l.strip()]
    rows = []
    t0 = time.time()
    for ri, d in enumerate(reps):
        nwk = os.path.join(d, f"tree_{TREEV}.nwk")
        if not os.path.exists(nwk):
            continue
        parent, elen, labels = read_newick_topology(nwk)
        st = {}
        with open(os.path.join(d, "states.csv")) as f:
            for r in csv.DictReader(f):
                st[r["cell"]] = (r["clone"], int(r["state"]), float(r["trait"]))
        scen = os.path.basename(d).split("_")[0]

        clones = sorted({v[0] for v in st.values()})
        for cl in clones:
            cells = [c for c, v in st.items() if v[0] == cl and c in labels]
            if len(cells) < 6:
                continue
            nodes = [labels[c] for c in cells]
            npar, nlen, leafpos = subtree_for(parent, elen, nodes)
            leaf_state = np.full(len(npar), -1, int)
            for c, nd in zip(cells, nodes):
                if nd in leafpos:
                    leaf_state[leafpos[nd]] = st[c][1]
            if TREEV == "unit":
                nlen = np.ones_like(nlen)
            tree = VarTree(npar, leaf_state, elen=nlen)
            t1 = time.time()
            try:
                fit = fit_var(tree)
            except Exception as e:
                rows.append(dict(rep_dir=d, scenario=scen, tree_variant=TREEV, clone=cl,
                                 n_cells=len(cells), StudyD_pi=np.nan, StudyD_k=np.nan,
                                 StudyD_k_at_bound="", StudyD_loglik=np.nan,
                                 StudyD_secs=time.time() - t1, StudyD_note=str(e)[:60]))
                continue
            rows.append(dict(rep_dir=d, scenario=scen, tree_variant=TREEV, clone=cl,
                             n_cells=len(cells), StudyD_pi=fit["pi"], StudyD_k=fit["k"],
                             StudyD_k_at_bound=int(fit["k_at_bound"]),
                             StudyD_loglik=fit["loglik"],
                             StudyD_secs=time.time() - t1, StudyD_note=""))
        if (ri + 1) % 5 == 0:
            print(f"  {ri+1}/{len(reps)} reps, {(time.time()-t0)/60:.1f} min", flush=True)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    cols = list(rows[0].keys())
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    fin = sum(1 for r in rows if np.isfinite(r["StudyD_k"]))
    print(f"wrote {OUT} with {len(rows)} rows; finite k: {fin}/{len(rows)}")


if __name__ == "__main__":
    main()
