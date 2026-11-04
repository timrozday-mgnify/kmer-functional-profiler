"""Tree index, phase 4's scoring (nearest.py): the ranking curve and earth mover's distance."""

import sys

import numpy as np
import pytest
from test_tree_phase0 import SCRIPT

sys.path.insert(0, str(SCRIPT.parent))
import nearest  # noqa: E402


def test_curve_counts_first_match_and_precision() -> None:
    # ranked units 0, 1, 2; unit 0 matches item 0, unit 2 items 0 and 1; item 2 is never found
    rec, prec = nearest.curve(
        np.array([0, 2, 2]), np.array([0, 0, 1]), np.arange(3, dtype=np.int64), 3
    )
    assert rec == pytest.approx([1 / 3, 1 / 3, 2 / 3])
    assert prec == pytest.approx([1, 1 / 2, 2 / 3])
    assert nearest.at_precision(rec, prec, 0.6) == pytest.approx(2 / 3)
    assert nearest.at_precision(rec, prec, 1.0) == pytest.approx(1 / 3)
    assert nearest.at_precision(rec, prec, 1.1) == 0


def test_emd() -> None:
    one = np.array([1.0])
    # the same sequence; a 90% neighbour; nothing in common
    i, j = np.array([0]), np.array([0])
    assert nearest.emd(one, one, i, j, np.array([1.0])) == pytest.approx(0)
    assert nearest.emd(one, one, i, j, np.array([0.9])) == pytest.approx(0.1)
    assert nearest.emd(one, one, i[:0], j[:0], np.zeros(0)) == 1
    # half the truth's mass is on a strain with no hit; mass is normalised first
    assert nearest.emd(np.array([2.0, 2.0]), np.array([3.0]), i, j, np.array([1.0])) == (
        pytest.approx(0.5)
    )
    # the LP takes the better pairing, not the greedy first one
    got = nearest.emd(
        np.array([1.0, 1.0]), np.array([1.0, 1.0]),
        np.array([0, 0, 1]), np.array([0, 1, 0]), np.array([0.9, 0.8, 0.8]),
    )  # fmt: skip
    assert got == pytest.approx(0.2)


def test_uses_sweep_and_restricted_forest() -> None:
    import uses

    # nodes over positions 0-7: (0,2), (2,4), (0,4), (4,8); only 1, 2, 3, 6 are kept
    left, right, tops, ranges = uses.restricted_forest(
        np.array([1, 2, 3, 6]), np.array([0, 2, 0, 4]), np.array([2, 4, 4, 8])
    )
    assert ranges[4:].tolist() == [[1, 3], [0, 3]]  # (0,2) and (4,8) hold one leaf each
    assert left.tolist() == [1, 0]
    assert right.tolist() == [2, 4]
    assert tops.tolist() == [5, 3]
    # k-mer 7: gained on leaves 0-3, lost on 1-2; k-mer 9: leaf 2
    h, leaf = uses.sweep(
        np.array([7, 7, 9], np.uint64), np.array([0, 1, 2]), np.array([4, 3, 3]),
        np.array([1, -1, 1], np.int8),
    )  # fmt: skip
    assert sorted(zip(h.tolist(), leaf.tolist(), strict=True)) == [(7, 0), (7, 3), (9, 2)]
