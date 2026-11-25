"""Benchmark IV: how much does tree reconstruction cost each method?

For each replicate we emit four trees over the same cells with the same states:

  tree_true.nwk    T1, the true sampled genealogy with true branch lengths
  tree_unit.nwk    T1 topology, unit edges
  tree_greedy.nwk  T3 from a greedy Cassiopeia-style reconstruction of simulated barcodes
  tree_nj.nwk      T3 from neighbour joining
  tree_upgma.nwk   T3 from UPGMA

Only the tree changes between them, so any difference in a method's output is attributable
to reconstruction, not to the data.

Usage: python -I reconstruction_bench.py <studyd_src> <outdir> <n_reps> <seed0> [recorder]
"""
import os, sys, json, csv
import numpy as np

STUDYD = sys.argv[1]
sys.path.insert(0, STUDYD)
from src.branching.recorder import (simulate_genealogy, induced_subtree, record_barcodes,
                                    reconstruct, saturation)

OUT, NREP, SEED0 = sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
RECORDER = sys.argv[5] if len(sys.argv) > 5 else "cassiopeia_like"

N_CLONES, PER_CLONE, N_GEN = 8, 48, 11
K_SHARED = 0.3
REC = dict(cassiopeia_like=dict(n_sites=30, mu=0.06, n_states=40, dropout=0.10, missing=0.10),
           gestalt_like=dict(n_sites=10, mu=0.20, n_states=12, dropout=0.15, missing=0.15))[RECORDER]


def evolve_binary(parent, rng, pi, k):
    E = np.exp(-k)
    st = np.empty(len(parent), int)
    st[0] = int(rng.random() < pi)
    for i in range(1, len(parent)):
        p = st[parent[i]]
        up = pi * (1 - E) if p == 0 else pi + (1 - pi) * E
        st[i] = int(rng.random() < up)
    return st


def evolve_ou(parent, rng, theta, alpha, sigma=1.0):
    x = np.empty(len(parent), float)
    x[0] = theta + rng.normal(0, sigma / np.sqrt(2 * alpha))
    e = np.exp(-alpha); v = sigma ** 2 * (1 - e ** 2) / (2 * alpha)
    for i in range(1, len(parent)):
        x[i] = theta + (x[parent[i]] - theta) * e + rng.normal(0, np.sqrt(v))
    return x


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


def graft(subtrees, stem=40.0):
    parent, elen, offsets = [-1], [0.0], []
    for sp, se in subtrees:
        off = len(parent); offsets.append(off)
        for i in range(len(sp)):
            parent.append(off + sp[i] if sp[i] >= 0 else 0)
            elen.append(float(se[i]) if sp[i] >= 0 else stem)
    return np.array(parent), np.array(elen), offsets


def main():
    os.makedirs(OUT, exist_ok=True)
    made, rep = 0, 0
    while made < NREP and rep < NREP * 6:
        rng = np.random.default_rng(SEED0 + rep); rep += 1
        per_clone, truth, sat_all, ok = [], [], [], True
        for c in range(N_CLONES):
            par, gen, live = simulate_genealogy(N_GEN, rng)
            if len(live) < PER_CLONE:
                ok = False; break
            take = rng.choice(live, PER_CLONE, replace=False)
            sub_par, sub_elen, leafmap, _ = induced_subtree(par, take)
            leaf_idx = sorted(leafmap.keys())
            orig = np.array([leafmap[i] for i in leaf_idx], int)

            pi_c = rng.uniform(0.15, 0.85)
            st = evolve_binary(par, rng, pi_c, K_SHARED)
            xt = evolve_ou(par, rng, 2.0 * (pi_c - 0.5), K_SHARED)

            # T2: record barcodes on the full genealogy, read them at the sampled cells
            A = record_barcodes(par, gen, rng, REC["n_sites"], REC["mu"], REC["n_states"],
                                REC["dropout"], REC["missing"])
            A_leaves = A[orig]
            sat_all.append(float(saturation(A_leaves)))
            recon = {}
            for meth in ("greedy", "nj", "upgma"):
                try:
                    rp = reconstruct(A_leaves, method=meth, rng=rng)
                except Exception:
                    rp = None
                recon[meth] = rp
            if any(v is None for v in recon.values()):
                ok = False; break
            per_clone.append(dict(par=sub_par, elen=sub_elen, leaf_idx=leaf_idx,
                                  state=st[orig], trait=xt[orig], clone=c, recon=recon,
                                  n_leaf=len(leaf_idx)))
            truth.append(dict(clone=c, pi=float(pi_c), k=K_SHARED,
                              theta=float(2.0 * (pi_c - 0.5)), alpha=K_SHARED))
        if not ok:
            continue

        d = os.path.join(OUT, f"R_rep{made:03d}")
        os.makedirs(d, exist_ok=True)

        # --- true sampled genealogy (T1)
        parent, elen, offsets = graft([(p["par"], p["elen"]) for p in per_clone])
        leaf_nodes, rows_csv = [], []
        for ct, off in zip(per_clone, offsets):
            for li, i in enumerate(ct["leaf_idx"]):
                leaf_nodes.append(off + i)
                rows_csv.append((f"cell_{len(leaf_nodes)-1}", f"C{ct['clone']}",
                                 int(ct["state"][li]), float(ct["trait"][li])))
        labels = [None] * len(parent)
        for j, v in enumerate(leaf_nodes):
            labels[v] = f"cell_{j}"
        open(os.path.join(d, "tree_true.nwk"), "w").write(to_newick(parent, elen, labels) + "\n")
        open(os.path.join(d, "tree_unit.nwk"), "w").write(
            to_newick(parent, np.ones(len(parent)), labels) + "\n")

        # --- reconstructed trees, leaves in the SAME order so labels line up
        for meth in ("greedy", "nj", "upgma"):
            subs = [(p["recon"][meth], np.ones(len(p["recon"][meth]))) for p in per_clone]
            rpar, relen, roff = graft(subs)
            rlabels = [None] * len(rpar)
            cursor = 0
            for ct, off in zip(per_clone, roff):
                sp = ct["recon"][meth]
                kids = set(int(x) for x in sp if x >= 0)
                leaves = [i for i in range(len(sp)) if i not in kids]
                for i in leaves:
                    rlabels[off + i] = f"cell_{cursor}"
                    cursor += 1
            if cursor != len(leaf_nodes):
                ok = False
                break
            open(os.path.join(d, f"tree_{meth}.nwk"), "w").write(
                to_newick(rpar, relen, rlabels) + "\n")
        if not ok:
            continue

        with open(os.path.join(d, "states.csv"), "w", newline="") as f:
            w = csv.writer(f); w.writerow(["cell", "clone", "state", "trait"])
            for r in rows_csv:
                w.writerow([r[0], r[1], r[2], f"{r[3]:.6f}"])
        json.dump(dict(recorder=RECORDER, rep=made, seed=SEED0 + rep - 1, clones=truth,
                       V_pi=float(np.var([t["pi"] for t in truth])), V_logk=0.0,
                       mean_k=K_SHARED, saturation=float(np.mean(sat_all))),
                  open(os.path.join(d, "truth.json"), "w"), indent=1)
        made += 1
    print(f"{RECORDER}: wrote {made} reconstruction replicates to {OUT}")


if __name__ == "__main__":
    main()
