"""Tests for the genealogy simulator and the observation chain T0 -> T1 -> T2 -> T3."""
import numpy as np
import pytest

from src.branching.recorder import (simulate_genealogy, induced_subtree, record_barcodes,
                                    reconstruct, saturation)


def test_genealogy_is_a_single_rooted_tree():
    rng = np.random.default_rng(0)
    par, gen, live = simulate_genealogy(8, rng)
    assert int((par < 0).sum()) == 1, "exactly one root"
    assert par[0] == -1
    assert len(live) > 0
    # every non-root node's parent precedes it, so the array is a valid topological order
    assert all(par[i] < i for i in range(1, len(par)))


def test_generation_depth_is_consistent():
    rng = np.random.default_rng(1)
    par, gen, live = simulate_genealogy(6, rng)
    for i in range(1, len(par)):
        assert gen[i] == gen[par[i]] + 1


def test_bottleneck_reduces_diversity():
    """A founder bottleneck should leave fewer distinct deep ancestors among live cells."""
    def n_anc(bott):
        rng = np.random.default_rng(4)
        par, gen, live = simulate_genealogy(9, rng, bottleneck_at=bott[0],
                                            bottleneck_size=bott[1]) if bott else \
              simulate_genealogy(9, rng)
        anc = set()
        for v in live:
            x = v
            while gen[x] > 3:
                x = par[x]
            anc.add(x)
        return len(anc)
    assert n_anc((3, 4)) <= n_anc(None)


def test_induced_subtree_preserves_sampled_cells():
    rng = np.random.default_rng(2)
    par, gen, live = simulate_genealogy(9, rng)
    take = rng.choice(live, 32, replace=False)
    sub_par, sub_elen, leafmap, _ = induced_subtree(par, take)
    assert int((sub_par < 0).sum()) == 1, "pruned tree keeps exactly one root"
    assert len(leafmap) == len(take), "every sampled cell survives as a leaf"
    assert set(leafmap.values()) == set(int(t) for t in take)
    assert (sub_elen >= 0).all()


def test_induced_subtree_edge_lengths_count_generations():
    """Collapsing unary chains must accumulate the generations they spanned."""
    rng = np.random.default_rng(6)
    par, gen, live = simulate_genealogy(8, rng)
    take = rng.choice(live, 16, replace=False)
    sub_par, sub_elen, leafmap, _ = induced_subtree(par, take)
    root = int(np.where(sub_par < 0)[0][0])
    for new_i, orig in leafmap.items():
        depth, x = 0.0, new_i
        while sub_par[x] >= 0:
            depth += sub_elen[x]
            x = sub_par[x]
        assert x == root
        assert depth == pytest.approx(gen[orig] - gen[int(np.where(par < 0)[0][0])], abs=1e-6)


def test_recorder_saturates_with_more_generations():
    rng = np.random.default_rng(3)
    sats = []
    for n_gen in (4, 10):
        par, gen, live = simulate_genealogy(n_gen, rng)
        A = record_barcodes(par, gen, rng, n_sites=10, mu=0.15, n_states=12,
                            dropout=0.0, missing=0.0)
        sats.append(float(saturation(A[live])))
    assert sats[1] >= sats[0], "a longer experiment edits more sites"


@pytest.mark.parametrize("method", ["greedy", "nj", "upgma"])
def test_reconstruction_returns_a_valid_tree(method):
    rng = np.random.default_rng(8)
    par, gen, live = simulate_genealogy(8, rng)
    take = rng.choice(live, 24, replace=False)
    A = record_barcodes(par, gen, rng, n_sites=20, mu=0.08, n_states=20,
                        dropout=0.05, missing=0.05)
    rp = reconstruct(A[take], method=method, rng=rng)
    rp = np.asarray(rp)
    assert int((rp < 0).sum()) == 1, f"{method} must return a single-rooted tree"
    kids = set(int(x) for x in rp if x >= 0)
    n_leaves = sum(1 for i in range(len(rp)) if i not in kids)
    assert n_leaves == 24, f"{method} must keep all 24 sampled cells as leaves"
