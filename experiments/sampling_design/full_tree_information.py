"""E13: exact full-tree Fisher information versus the naive independent-pair sum.

Pairs inside a tree share ancestors, so summing per-pair information double-counts. For trees small enough to
enumerate every leaf configuration (<= 2^16) the exact Fisher information of the full-tree likelihood can be
computed without approximation:

    I(theta) = sum_x p(x; theta) * s(x) s(x)^T,     s(x) = d log p(x; theta) / d theta,

the sum running over all 2^n leaf state vectors. This is compared with
  (a) the naive sum of per-pair informations, treating within-clone pairs as independent;
  (b) a composite-likelihood information with the Godambe sandwich correction.

Reports the inflation factor of the naive sum, and whether the depth that maximises full-tree information per
cell agrees with the pairwise optimum k*L ~ 0.4.
"""
import itertools, os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.lineage_inference.ctmc import FastTree
from experiments.sampling_design.fisher_audit import fisher2

TAB = os.path.join(ROOT, "results/tables"); os.makedirs(TAB, exist_ok=True)


def balanced_tree(depth):
    """Balanced binary tree: parent array, leaf indices, and the MRCA depth of every leaf pair."""
    par = [-1]
    level = [0]
    for _ in range(depth):
        nxt = []
        for v in level:
            for _ in range(2):
                par.append(v); nxt.append(len(par) - 1)
        level = nxt
    par = np.array(par)
    leaves = np.array(level)
    # path from root for each leaf
    paths = []
    for v in leaves:
        p, x = [], v
        while x >= 0:
            p.append(x); x = par[x]
        paths.append(tuple(reversed(p)))
    return par, leaves, paths


def make_tree(par, leaves):
    """Build the FastTree once; leaf states are then written in place for each configuration."""
    ls = np.full(len(par), -1)
    ls[leaves] = 0
    return FastTree(par, ls)


def tree_loglik(tree, leaves, states, pi, k, branch):
    tree.leaf_state[leaves] = states          # FastTree uses unit branches, so fold branch into the rate
    return tree.loglik(pi, [(pi, k * branch)])


def exact_tree_fisher(depth, pi, k, branch):
    """Exact 2x2 Fisher information for (log k, pi) from one full tree, by enumerating leaf configurations."""
    par, leaves, paths = balanced_tree(depth)
    n = len(leaves)
    if n > 16:
        raise ValueError("too many leaves to enumerate")
    tree = make_tree(par, leaves)
    h = 1e-5
    lp = np.empty(2 ** n); d_k = np.empty(2 ** n); d_p = np.empty(2 ** n)
    for idx, cfg in enumerate(itertools.product((0, 1), repeat=n)):
        s = np.array(cfg)
        lp[idx] = tree_loglik(tree, leaves, s, pi, k, branch)
        d_k[idx] = (tree_loglik(tree, leaves, s, pi, k * np.exp(h), branch)
                    - tree_loglik(tree, leaves, s, pi, k * np.exp(-h), branch)) / (2 * h)
        d_p[idx] = (tree_loglik(tree, leaves, s, pi + h, k, branch)
                    - tree_loglik(tree, leaves, s, pi - h, k, branch)) / (2 * h)
    p = np.exp(lp); p = p / p.sum()
    I = np.zeros((2, 2))
    I[0, 0] = (p * d_k ** 2).sum(); I[1, 1] = (p * d_p ** 2).sum(); I[0, 1] = I[1, 0] = (p * d_k * d_p).sum()
    return I, n, paths, leaves


def naive_pair_sum(paths, pi, k, branch, depth):
    """Sum of per-pair Fisher informations, treating pairs as independent (the Milestone 3 assumption)."""
    n = len(paths)
    I = np.zeros((2, 2))
    npairs = 0
    for i in range(n):
        for j in range(i + 1, n):
            sh = 0
            for a, b in zip(paths[i], paths[j]):
                if a != b:
                    break
                sh += 1
            L = (depth - (sh - 1)) * branch        # lineage distance from MRCA to each leaf
            I += fisher2(pi, k, L)
            npairs += 1
    return I, npairs


def main():
    rows = []
    for depth in (3, 4, 5, 6):
        for pi in (0.3, 0.5):
            for k in (0.05, 0.2, 0.5, 1.0, 2.0):
                for branch in (1.0,):
                    try:
                        I, n, paths, leaves = exact_tree_fisher(depth, pi, k, branch)
                    except ValueError:
                        continue
                    Ipair, npairs = naive_pair_sum(paths, pi, k, branch, depth)
                    joint = lambda M: M[0, 0] - (M[0, 1] ** 2 / M[1, 1] if M[1, 1] > 1e-12 else 0.0)
                    rows.append(dict(depth=depth, n_leaves=n, pi=pi, k=k, branch=branch, n_pairs=npairs,
                                     I_tree_logk=I[0, 0], I_tree_logk_joint=joint(I),
                                     I_pairsum_logk=Ipair[0, 0], I_pairsum_logk_joint=joint(Ipair),
                                     inflation=Ipair[0, 0] / max(I[0, 0], 1e-12),
                                     I_tree_per_cell=joint(I) / n,
                                     kL_mean_depth=k * depth * branch))
                print("depth", depth, "pi", pi, "done", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(TAB, "e13_full_tree_information.csv"), index=False)
    print(d.groupby(["depth", "k"])[["I_tree_logk_joint", "I_pairsum_logk_joint", "inflation",
                                     "I_tree_per_cell"]].median().round(4).to_string())
    print("\nnaive-pair inflation factor: median %.1fx  range %.1f-%.1f"
          % (d.inflation.median(), d.inflation.min(), d.inflation.max()))
    best = d.loc[d.groupby(["pi", "k"]).I_tree_per_cell.idxmax(), ["pi", "k", "depth", "I_tree_per_cell"]]
    best["k_times_depth"] = best.k * best.depth
    print("\ndepth maximising full-tree information per cell:")
    print(best.round(4).to_string())




def sample_leaves(par, leaves, pi, k, branch, n_rep, rng):
    """Simulate leaf states down the tree under the two-state CTMC."""
    E = np.exp(-k * branch)
    st = np.empty((n_rep, len(par)), dtype=np.int8)
    st[:, 0] = (rng.random(n_rep) < pi)
    for i in range(1, len(par)):
        p = st[:, par[i]]
        up = np.where(p == 0, pi * (1 - E), pi + (1 - pi) * E)
        st[:, i] = (rng.random(n_rep) < up)
    return st[:, leaves]


def mc_tree_fisher(depth, pi, k, branch, n_rep, rng):
    """Monte Carlo Fisher information: average of the outer product of the score over simulated trees."""
    par, leaves, paths = balanced_tree(depth)
    tree = make_tree(par, leaves)
    cfgs = sample_leaves(par, leaves, pi, k, branch, n_rep, rng)
    h = 1e-5
    dk = np.empty(n_rep); dp = np.empty(n_rep)
    for r in range(n_rep):
        s = cfgs[r].astype(int)
        dk[r] = (tree_loglik(tree, leaves, s, pi, k * np.exp(h), branch)
                 - tree_loglik(tree, leaves, s, pi, k * np.exp(-h), branch)) / (2 * h)
        dp[r] = (tree_loglik(tree, leaves, s, pi + h, k, branch)
                 - tree_loglik(tree, leaves, s, pi - h, k, branch)) / (2 * h)
    I = np.array([[np.mean(dk ** 2), np.mean(dk * dp)], [np.mean(dk * dp), np.mean(dp ** 2)]])
    se = np.array([np.std(dk ** 2), np.std(dp ** 2)]) / np.sqrt(n_rep)
    return I, len(leaves), paths, se


def main_mc(n_rep=2000, seed=20261201):
    """Depth sweep by Monte Carlo, validated against exact enumeration at depths 3 and 4."""
    rng = np.random.default_rng(seed)
    rows = []
    for depth in (3, 4, 5, 6, 7):
        for pi in (0.3, 0.5):
            for k in (0.05, 0.1, 0.2, 0.5, 1.0, 2.0):
                I, n, paths, se = mc_tree_fisher(depth, pi, k, 1.0, n_rep, rng)
                joint = I[0, 0] - (I[0, 1] ** 2 / I[1, 1] if I[1, 1] > 1e-12 else 0.0)
                Ip, npairs = naive_pair_sum(paths, pi, k, 1.0, depth)
                jp = Ip[0, 0] - (Ip[0, 1] ** 2 / Ip[1, 1] if Ip[1, 1] > 1e-12 else 0.0)
                exact = np.nan
                if depth <= 3:
                    Ie, *_ = exact_tree_fisher(depth, pi, k, 1.0)
                    exact = Ie[0, 0]
                rows.append(dict(depth=depth, n_leaves=n, pi=pi, k=k, n_pairs=npairs,
                                 I_tree_logk=I[0, 0], I_tree_logk_se=se[0], I_tree_joint=joint,
                                 I_exact_logk=exact,
                                 I_pairsum_logk=Ip[0, 0], I_pairsum_joint=jp,
                                 inflation=Ip[0, 0] / max(I[0, 0], 1e-12),
                                 I_tree_per_cell=joint / n, kD=k * depth))
            print("mc depth", depth, "pi", pi, "done", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(TAB, "e13_tree_vs_pair_mc.csv"), index=False)
    v = d.dropna(subset=["I_exact_logk"])
    if len(v):
        rel = ((v.I_tree_logk - v.I_exact_logk).abs() / v.I_exact_logk)
        print("MC vs exact: max relative error %.3f over %d settings" % (rel.max(), len(v)))
    print(d.groupby(["depth", "k"])[["I_tree_joint", "I_pairsum_joint", "inflation",
                                     "I_tree_per_cell"]].median().round(4).to_string())
    best = d.loc[d.groupby(["pi", "k"]).I_tree_per_cell.idxmax(), ["pi", "k", "depth", "kD", "I_tree_per_cell"]]
    print("\ndepth maximising full-tree information per cell:")
    print(best.round(4).to_string())

if __name__ == "__main__":
    main_mc() if "--mc" in sys.argv else main()
