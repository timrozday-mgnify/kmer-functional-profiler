"""Tree index, tier 2 (tier2.py): chunks along the leaf order and the routed work list."""

import sys

import numpy as np
import polars as pl
from test_tree_phase0 import SCRIPT

sys.path.insert(0, str(SCRIPT.parent))
import tier2  # noqa: E402


def test_chunks_pack_leaves_in_order_under_the_budget() -> None:
    length = np.array([5, 5, 5, 12, 3, 3])
    leaves = np.array([0, 1, 2, 3, 4, 5])
    # 5+5 fit 10; 5 starts a new chunk; 12 is alone (longer than the budget); 3+3
    assert tier2.chunks(leaves, length, 10).tolist() == [0, 0, 1, 2, 3, 3]


def test_work_list_routes_reads_to_groups_above_the_floor() -> None:
    group = np.array([0, 0, 2, 2, 4])
    aligned = pl.DataFrame(
        {"read": [0, 0, 1, 1], "pos": [1, 3, 4, 0], "score": [30.0, 12.0, 25.0, 5.0]}
    )
    got = tier2.work_list(aligned, group, 20.0).sort("read", "g")
    assert got.rows() == [(0, 0), (1, 4)]  # read 0's pair at 12 bits and read 1's at 5 dropped


def test_em_gives_shared_reads_to_the_supported_leaf() -> None:
    """Reads tied between leaves 1 and 2 go to 1, which other reads support alone."""
    calls = pl.DataFrame(
        {"read": [0, 1, 2, 3, 3, 4, 4], "pos": [1, 1, 1, 1, 2, 1, 2], "score": [30.0] * 7}
    )
    w = tier2.em_weights(calls, 30)
    assert abs(w.group_by("read").agg(pl.col("weight").sum())["weight"] - 1).max() < 1e-9
    shared = w.filter(pl.col("read") == 3).sort("pos")["weight"].to_list()
    assert shared[0] > 0.95 > 0.05 > shared[1]


def test_budget_splits_memory() -> None:
    residues, pairs = tier2.budget(3.0)
    assert residues * tier2.RESIDUE_BYTES == pairs * tier2.PAIR_BYTES == 2**30


def test_read_scoring_pools_small_depth_bins() -> None:
    import recall

    # 10+40 reach 30; 5+30 do; the short tail (3) joins the last full pool
    pools = recall.merge_bins([10, 40, 5, 30, 3], 30)
    assert pools == [0, 0, 1, 1, 1]
    assert [recall.span(pools, p) for p in (0, 1)] == ["0-1", "1-"]
    assert recall.merge_bins([5, 5], 30) == [0, 0]  # too few anywhere: one pool
    m, lo, hi = recall.mean_ci(np.array([0.0, 1.0] * 50))
    assert m == 0.5 and lo < 0.5 < hi


def test_route_k_grows_with_the_candidates() -> None:
    assert [tier2.route_k(int(r)) for r in (2.5e7, 6e8, 1e12)] == [5, 6, 7]
