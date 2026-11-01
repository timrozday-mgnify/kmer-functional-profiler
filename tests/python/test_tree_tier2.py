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
