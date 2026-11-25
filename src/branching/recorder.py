"""Recorder-realistic generative benchmark: the T0 -> T1 -> T2 -> T3 observation hierarchy.

T0  true genealogy, every cell including those never sampled
T1  sampled genealogy, the induced subtree on sampled cells (ancestry known without error)
T2  recorded information: the allele matrix the recorder actually writes, with saturation, homoplasy,
    dropout and missing data
T3  reconstructed genealogy from T2 by a named algorithm

Reconstruction algorithms implemented here are standard and are stand-ins for the published pipelines, which
are not installable in this environment:
  * `greedy`  - character-based recursive splitting on the most frequent non-ancestral state. This is the
                algorithm Cassiopeia's greedy solver uses, reimplemented from its description.
  * `nj`      - neighbour joining on normalised Hamming distance over shared non-missing sites.
  * `upgma`   - average-linkage on the same distance.
Parameter provenance for the two recorder settings is in reports/milestone5_recorder_simulation.md.
"""
import numpy as np

MISSING = -1
UNEDITED = 0


# --------------------------------------------------------------------------------------- T0 and T1

def simulate_genealogy(n_gen, rng, p_death=0.0, clone_bias=None, bottleneck_at=None, bottleneck_size=None):
    """Synchronous-division genealogy for n_gen generations.

    clone_bias: optional per-founder relative division probability, producing unequal clonal expansion.
    bottleneck_at/size: at that generation, keep only `bottleneck_size` randomly chosen lineages, modelling a
    metastatic founder bottleneck; this is the main source of genealogical imbalance.
    Returns parent array, generation of each node, and the indices of the final live cells.
    """
    parent = [-1]
    gen = [0]
    live = [0]
    for g in range(1, n_gen + 1):
        if bottleneck_at is not None and g == bottleneck_at and bottleneck_size < len(live):
            live = list(rng.choice(live, bottleneck_size, replace=False))
        nxt = []
        for v in live:
            pdiv = 1.0 if clone_bias is None else clone_bias(v)
            if rng.random() > pdiv:                       # quiescent this generation
                nxt.append(v); continue
            for _ in range(2):
                if rng.random() < p_death:
                    continue
                parent.append(v); gen.append(g); nxt.append(len(parent) - 1)
        live = nxt
        if not live:
            break
    return np.array(parent), np.array(gen), np.array(live)


def simulate_genealogy_tracked(n_gen, rng, p_death=0.0, bottleneck_at=None, bottleneck_size=None):
    """As simulate_genealogy, but also returns the set of live cells at the END of every generation.

    Needed for multi-timepoint sampling designs, where cells are harvested at more than one time.
    """
    parent = [-1]; gen = [0]; live = [0]; snapshots = [list(live)]
    for g in range(1, n_gen + 1):
        if bottleneck_at is not None and g == bottleneck_at and bottleneck_size < len(live):
            live = list(rng.choice(live, bottleneck_size, replace=False))
        nxt = []
        for v in live:
            for _ in range(2):
                if rng.random() < p_death:
                    continue
                parent.append(v); gen.append(g); nxt.append(len(parent) - 1)
        live = nxt
        snapshots.append(list(live))
        if not live:
            break
    return np.array(parent), np.array(gen), np.array(live), snapshots


def evolve_states(parent, rng, pi, k):
    """Two-state CTMC down the genealogy, one unit of time per generation step."""
    E = np.exp(-k)
    st = np.empty(len(parent), int)
    st[0] = int(rng.random() < pi)
    for i in range(1, len(parent)):
        p = st[parent[i]]
        up = pi * (1 - E) if p == 0 else pi + (1 - pi) * E
        st[i] = int(rng.random() < up)
    return st


def induced_subtree(parent, keep):
    """T1: the genealogy induced on `keep`, with unary chains collapsed.

    Returns a new parent array, the mapping from new leaf index to original cell, and the number of original
    generations spanned by each new edge (true branch length in divisions).
    """
    keep = list(keep)
    needed = set(keep)
    for v in keep:
        x = v
        while parent[x] >= 0:
            x = parent[x]
            needed.add(x)
    # a node is retained if it is a kept leaf or has >= 2 retained children
    children = {}
    for v in sorted(needed):
        p = parent[v]
        if p >= 0 and p in needed:
            children.setdefault(p, []).append(v)
    keepset = set(keep)

    def n_branching(v):
        return len(children.get(v, []))

    retained = {v for v in needed if v in keepset or n_branching(v) >= 2}
    root = min(needed)
    retained.add(root)
    newpar, length, order = {}, {}, sorted(retained)
    idx = {v: i for i, v in enumerate(order)}
    for v in order:
        if v == root:
            newpar[v] = -1; length[v] = 0; continue
        x = parent[v]; L = 1
        while x not in retained:
            x = parent[x]; L += 1
        newpar[v] = idx[x]; length[v] = L
    par = np.array([newpar[v] for v in order])
    ln = np.array([length[v] for v in order], float)
    leafmap = {idx[v]: v for v in keep if v in retained}
    return par, ln, leafmap, order


# --------------------------------------------------------------------------------------- T2

def record_barcodes(parent, gen, rng, n_sites, mu, n_states, dropout, missing,
                    state_weights=None):
    """Cumulative recorder: each unedited site edits with probability mu per generation.

    n_states    number of distinct indel outcomes per site (homoplasy: the same outcome recurs independently)
    state_weights  probability over outcomes; a skewed distribution raises homoplasy
    dropout     probability a site is silenced for a whole lineage once it occurs (modelled per cell as loss)
    missing     probability a site is unobserved in the readout of a given cell
    Returns the allele matrix for every node (cells x sites).
    """
    if state_weights is None:
        w = 1.0 / np.arange(1, n_states + 1) ** 1.0      # Zipf-like outcome frequencies
        state_weights = w / w.sum()
    A = np.zeros((len(parent), n_sites), int)
    for i in range(1, len(parent)):
        a = A[parent[i]].copy()
        un = np.where(a == UNEDITED)[0]
        if len(un):
            hit = un[rng.random(len(un)) < mu]
            if len(hit):
                a[hit] = rng.choice(np.arange(1, n_states + 1), len(hit), p=state_weights)
        if dropout > 0:
            d = rng.random(n_sites) < dropout
            a[d & (a != MISSING)] = MISSING
        A[i] = a
    if missing > 0:
        M = rng.random(A.shape) < missing
        A = A.copy(); A[M] = MISSING
    return A


def saturation(A_leaves):
    """Fraction of observed sites that are already edited (1.0 means no further information can accrue)."""
    obs = A_leaves != MISSING
    if obs.sum() == 0:
        return np.nan
    return float((A_leaves[obs] > 0).mean())


# --------------------------------------------------------------------------------------- T3

def _dist_matrix(A):
    n = len(A)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            m = (A[i] != MISSING) & (A[j] != MISSING)
            D[i, j] = D[j, i] = (A[i][m] != A[j][m]).mean() if m.sum() else 1.0
    return D


def reconstruct(A, method="greedy", rng=None):
    """Return a parent array over 2n-1 nodes with the n leaves first, from the allele matrix."""
    n = len(A)
    if method in ("nj", "upgma"):
        return _distance_tree(_dist_matrix(A), method)
    return _greedy_tree(A, rng or np.random.default_rng(0))


def _distance_tree(D, method):
    n = len(D)
    active = list(range(n))
    par = [-1] * n
    D = D.copy().astype(float)
    sizes = {i: 1 for i in range(n)}
    while len(active) > 1:
        best, bi, bj = np.inf, None, None
        if method == "nj" and len(active) > 2:
            r = {i: sum(D[i, j] for j in active if j != i) / (len(active) - 2) for i in active}
            for ii, i in enumerate(active):
                for j in active[ii + 1:]:
                    q = D[i, j] - r[i] - r[j]
                    if q < best:
                        best, bi, bj = q, i, j
        else:
            for ii, i in enumerate(active):
                for j in active[ii + 1:]:
                    if D[i, j] < best:
                        best, bi, bj = D[i, j], i, j
        new = len(par); par.append(-1)
        par[bi] = new; par[bj] = new
        nd = np.zeros((new + 1, new + 1))
        nd[:D.shape[0], :D.shape[1]] = D
        for m in active:
            if m in (bi, bj):
                continue
            if method == "upgma":
                v = (D[bi, m] * sizes[bi] + D[bj, m] * sizes[bj]) / (sizes[bi] + sizes[bj])
            else:
                v = 0.5 * (D[bi, m] + D[bj, m] - D[bi, bj])
            nd[new, m] = nd[m, new] = v
        D = nd
        sizes[new] = sizes[bi] + sizes[bj]
        active = [m for m in active if m not in (bi, bj)] + [new]
    return np.array(par)


def _greedy_tree(A, rng):
    """Cassiopeia-style greedy: recursively split on the most frequent non-ancestral character state."""
    n = len(A)
    par = [-1] * n
    nodes = {}

    def split(idx, parent_node):
        if len(idx) == 1:
            par[idx[0]] = parent_node
            return
        best, bc, bs = -1, None, None
        for c in range(A.shape[1]):
            col = A[idx, c]
            for s in np.unique(col[col > 0]):
                cnt = int((col == s).sum())
                if 0 < cnt < len(idx) and cnt > best:
                    best, bc, bs = cnt, c, s
        if bc is None:                       # no informative character: star polytomy, resolved arbitrarily
            left, right = idx[:len(idx) // 2], idx[len(idx) // 2:]
        else:
            m = A[idx, bc] == bs
            left = [i for i, t in zip(idx, m) if t]; right = [i for i, t in zip(idx, m) if not t]
        new = len(par); par.append(parent_node)
        nodes[new] = True
        for part in (left, right):
            if len(part) == 1:
                par[part[0]] = new
            else:
                split(list(part), new)

    root = len(par); par.append(-1)
    split(list(range(n)), root)
    return np.array(par)


# --------------------------------------------------------------------------------------- metrics

def mrca_depths(par, n_leaves, max_pairs=20000, rng=None):
    """Shared-ancestor depth for leaf pairs, where depth counts shared nodes from the root."""
    paths = []
    for v in range(n_leaves):
        p, x = [], v
        while x >= 0:
            p.append(x); x = par[x]
        paths.append(tuple(reversed(p)))
    rng = rng or np.random.default_rng(0)
    pairs = [(i, j) for i in range(n_leaves) for j in range(i + 1, n_leaves)]
    if len(pairs) > max_pairs:
        pairs = [pairs[t] for t in rng.choice(len(pairs), max_pairs, replace=False)]
    out = []
    for i, j in pairs:
        a, b = paths[i], paths[j]
        s = 0
        for x, y in zip(a, b):
            if x != y:
                break
            s += 1
        out.append(s)
    return np.array(out), paths


def colless(par, n_leaves):
    """Colless imbalance, normalised; higher means more star-like / caterpillar-like."""
    children = {}
    for v, p in enumerate(par):
        if p >= 0:
            children.setdefault(p, []).append(v)
    nleaf = {}

    def count(v):
        if v not in children:
            return 1
        c = sum(count(u) for u in children[v])
        nleaf[v] = c
        return c

    root = int(np.where(np.asarray(par) < 0)[0][0])
    count(root)
    tot = 0
    for v, ch in children.items():
        if len(ch) == 2:
            a = nleaf.get(ch[0], 1); b = nleaf.get(ch[1], 1)
            tot += abs(a - b)
    denom = (n_leaves - 1) * (n_leaves - 2) / 2
    return tot / denom if denom > 0 else np.nan
