"""Benchmark II: vary the distribution of ancestral distances and see which methods track it.

Everything except tree geometry is held fixed: 12 clones, occupancy drawn from the same
distribution, and a single shared switching rate k = 0.3 per generation. Any apparent
rate heterogeneity a method reports is therefore an artefact of the geometry.

Geometries
  G_shallow    7 generations, random sampling
  G_deep      13 generations, random sampling
  G_clustered 11 generations, cells drawn from a few whole clades (close relatives common)
  G_dispersed 11 generations, cells spread across distinct clades (close relatives scarce)
  G_starlike  11 generations with an early bottleneck, which flattens the genealogy

Usage: python -I geometry_bench.py <studyd_src> <outdir> <geometry> <n_reps> <seed0>
"""
import os, sys, json, csv
import numpy as np

STUDYD = sys.argv[1]
sys.path.insert(0, STUDYD)
from src.branching.recorder import simulate_genealogy, induced_subtree

OUT, GEO, NREP, SEED0 = sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])

N_CLONES = 12
PER_CLONE = 48
K_SHARED = 0.3

CFG = {
    "G_shallow":   dict(n_gen=7,  sampling="random",    bott=None),
    "G_deep":      dict(n_gen=13, sampling="random",    bott=None),
    "G_clustered": dict(n_gen=11, sampling="clustered", bott=None),
    "G_dispersed": dict(n_gen=11, sampling="dispersed", bott=None),
    "G_starlike":  dict(n_gen=11, sampling="random",    bott=(3, 4)),
}
if GEO not in CFG:
    raise SystemExit(f"unknown geometry {GEO}; choose from {list(CFG)}")
cfg = CFG[GEO]


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
    e = np.exp(-alpha)
    v = (sigma ** 2) * (1 - e ** 2) / (2 * alpha)
    for i in range(1, len(parent)):
        x[i] = theta + (x[parent[i]] - theta) * e + rng.normal(0, np.sqrt(v))
    return x


def pick_cells(par, gen, live, rng, how, n):
    """Choose which live cells to sample, which is what sets the ancestral-distance spread."""
    live = np.asarray(live)
    if len(live) <= n:
        return live
    if how == "random":
        return rng.choice(live, n, replace=False)
    # group live cells by their ancestor at an intermediate generation
    target = max(1, gen.max() // 2)
    anc = {}
    for v in live:
        x = v
        while gen[x] > target:
            x = par[x]
        anc.setdefault(x, []).append(v)
    groups = list(anc.values())
    rng.shuffle(groups)
    if how == "clustered":
        out = []
        for g in groups:
            out.extend(g)
            if len(out) >= n:
                break
        return np.array(out[:n])
    if how == "dispersed":
        out, i = [], 0
        while len(out) < n:
            progressed = False
            for g in groups:
                if i < len(g):
                    out.append(g[i]); progressed = True
                    if len(out) >= n:
                        break
            if not progressed:
                break
            i += 1
        return np.array(out[:n])
    raise ValueError(how)


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


def mrca_depth_stats(parent, elen, leaf_idx):
    """Realised ancestral-distance distribution over sampled leaf pairs."""
    anc = {}
    for v in leaf_idx:
        chain, x, acc = {}, v, 0.0
        while x >= 0:
            chain[x] = acc
            if parent[x] < 0:
                break
            acc += elen[x]
            x = parent[x]
        anc[v] = chain
    ds = []
    L = list(leaf_idx)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            a, b = anc[L[i]], anc[L[j]]
            common = set(a) & set(b)
            if not common:
                continue
            ds.append(min(max(a[c], b[c]) for c in common))
    return np.array(ds) if ds else np.array([np.nan])


def main():
    os.makedirs(OUT, exist_ok=True)
    made, rep = 0, 0
    while made < NREP and rep < NREP * 6:
        rng = np.random.default_rng(SEED0 + rep); rep += 1
        clone_trees, truth = [], []
        ok = True
        for c in range(N_CLONES):
            ba, bs = (cfg["bott"] if cfg["bott"] else (None, None))
            par, gen, live = simulate_genealogy(cfg["n_gen"], rng, bottleneck_at=ba,
                                                bottleneck_size=bs)
            if len(live) < 8:
                ok = False; break
            take = pick_cells(par, gen, live, rng, cfg["sampling"], PER_CLONE)
            sub_par, sub_elen, leafmap, _ = induced_subtree(par, take)
            leaf_idx = sorted(leafmap.keys())
            orig = np.array([leafmap[i] for i in leaf_idx], int)
            pi_c = rng.uniform(0.15, 0.85)
            st = evolve_binary(par, rng, pi_c, K_SHARED)
            xt = evolve_ou(par, rng, 2.0 * (pi_c - 0.5), K_SHARED)
            clone_trees.append(dict(parent=sub_par, elen=sub_elen, leaf_idx=leaf_idx,
                                    state=st[orig], trait=xt[orig], clone=c))
            truth.append(dict(clone=c, pi=float(pi_c), k=K_SHARED,
                              theta=float(2.0 * (pi_c - 0.5)), alpha=K_SHARED))
        if not ok:
            continue

        parent, elen, offsets = [-1], [0.0], []
        for ct in clone_trees:
            off = len(parent); offsets.append(off)
            sp, se = ct["parent"], ct["elen"]
            for i in range(len(sp)):
                parent.append(off + sp[i] if sp[i] >= 0 else 0)
                elen.append(float(se[i]) if sp[i] >= 0 else 40.0)
        parent, elen = np.array(parent), np.array(elen)

        leaf_nodes, rows_csv = [], []
        for ct, off in zip(clone_trees, offsets):
            for li, i in enumerate(ct["leaf_idx"]):
                leaf_nodes.append(off + i)
                rows_csv.append((f"cell_{len(leaf_nodes)-1}", f"C{ct['clone']}",
                                 int(ct["state"][li]), float(ct["trait"][li])))
        labels = [None] * len(parent)
        for j, v in enumerate(leaf_nodes):
            labels[v] = f"cell_{j}"

        d = os.path.join(OUT, f"{GEO}_rep{made:03d}")
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, "tree_true.nwk"), "w").write(to_newick(parent, elen, labels) + "\n")
        open(os.path.join(d, "tree_unit.nwk"), "w").write(
            to_newick(parent, np.ones(len(parent)), labels) + "\n")
        with open(os.path.join(d, "states.csv"), "w", newline="") as f:
            w = csv.writer(f); w.writerow(["cell", "clone", "state", "trait"])
            for r in rows_csv:
                w.writerow([r[0], r[1], r[2], f"{r[3]:.6f}"])

        # realised geometry, measured within clones (the pairs the methods actually use)
        dd = []
        for ct, off in zip(clone_trees, offsets):
            idx = [off + i for i in ct["leaf_idx"]]
            dd.append(mrca_depth_stats(parent, elen, idx))
        dd = np.concatenate(dd)
        dd = dd[np.isfinite(dd)]
        json.dump(dict(geometry=GEO, rep=made, seed=SEED0 + rep - 1,
                       n_gen=cfg["n_gen"], sampling=cfg["sampling"],
                       clones=truth, V_pi=float(np.var([t["pi"] for t in truth])),
                       V_logk=0.0, mean_k=K_SHARED,
                       mrca_median=float(np.median(dd)), mrca_q25=float(np.percentile(dd, 25)),
                       mrca_q75=float(np.percentile(dd, 75)),
                       frac_close=float(np.mean(dd <= 2.0)),
                       kL_median=float(K_SHARED * np.median(dd)),
                       frac_in_window=float(np.mean((K_SHARED * dd >= 0.10) &
                                                    (K_SHARED * dd <= 1.0)))),
                  open(os.path.join(d, "truth.json"), "w"), indent=1)
        made += 1
    print(f"{GEO}: wrote {made} replicates")


if __name__ == "__main__":
    main()
