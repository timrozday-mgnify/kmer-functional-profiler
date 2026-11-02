"""Tree index: the leaf-ordered store (store.py) reads what a members scan would."""

import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from test_tree_phase0 import SCRIPT

sys.path.insert(0, str(SCRIPT.parent))
import store  # noqa: E402


def members_and_parts(tmp: Path) -> tuple[list[str], list[tuple[Path, int]]]:
    """Two parts' proteins (part-local positions, one with two proteins and a leafless -1),
    and their members."""
    m = pl.DataFrame(
        {
            "protein_id": [1, 2, 3, 4, 5, 6],
            "cluster_rep": [1, 1, 3, 4, 4, 6],
            "full_length": [True] * 6,
            "sequence": ["MKV", "MKVLA", "QQ", "WWWW", "HH", "PP"],
        }
    )
    m.write_parquet(tmp / "members.parquet")
    parts = []
    for name, ids, pos, offset in (("a", [1, 2, 3], [0, 0, 1], 0), ("b", [4, 5, 6], [1, 0, -1], 2)):
        f = tmp / f"{name}.proteins.parquet"
        pl.DataFrame({"protein_id": ids, "cluster_rep": [1, 1, 3] if name == "a" else [4, 4, 6],
                      "pos": pos}).write_parquet(f)  # fmt: skip
        parts.append((f, offset))
    return [str(tmp / "members.parquet")], parts


@pytest.mark.parametrize("fmt", store.FORMATS)
def test_store_reads_each_positions_longest_protein(tmp_path: Path, fmt: str) -> None:
    members, parts = members_and_parts(tmp_path)
    out = tmp_path / "store"
    for (f, offset), name in zip(parts, ("a", "b"), strict=True):
        store.build(f, members, offset, out, name, fmt, row_group=1)
    ls = store.LeafStore(out)
    got = ls.get(np.array([3, 0, 2, 1, 9]), ids=True)  # 9: no such position
    assert got["pos"].to_list() == [0, 1, 2, 3]
    # position 0 holds proteins 1 and 2: the longer; part b's positions start at 2
    assert got["sequence"].to_list() == ["MKVLA", "QQ", "HH", "WWWW"]
    assert got["protein_id"].to_list() == [2, 3, 5, 4]
    assert got["cluster_rep"].to_list() == [1, 3, 4, 4]
