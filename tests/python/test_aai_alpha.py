"""The aai-alpha workflow's survival counting and fit."""

import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from kmer_functional_profiler.query import survival

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflows" / "aai-alpha"))
import alpha  # noqa: E402

K = 3


def test_pair_survival_leaves_the_query_out_of_its_own_cluster() -> None:
    # Cluster 1: members 1 (ABCDEFG), 2 (ABCDEFH); cluster 5: member 5 (ABCDXYZ). P = 2 in its
    # own cluster 1 and against cluster 5; P = 9 (a fragment-free outsider) against 1.
    members = pl.DataFrame(
        {
            "protein_id": [1, 2, 5, 9],
            "cluster_rep": [1, 1, 5, 9],
            "full_length": [True, True, True, True],
            "sequence": ["MKVLAWE", "MKVLAWQ", "MKVLRST", "MKVLAWE"],
        }
    )
    pairs = pl.DataFrame(
        {"protein_id": [2, 2, 9], "cluster_rep": [1, 5, 1], "identity": [0.86, 0.57, 1.0]}
    )
    got = alpha.pair_survival(pairs, members, K, "protein").sort("protein_id", "cluster_rep")
    rows = {(r["protein_id"], r["cluster_rep"]): r for r in got.iter_rows(named=True)}
    # 2 vs 1 without itself: member 1 alone, 5 k-mers; 2 shares MKV KVL VLA LAW (4)
    assert rows[2, 1]["shared"] == 4 and rows[2, 1]["pin_sum"] == 5 and rows[2, 1]["n_members"] == 1
    assert rows[2, 5]["shared"] == 2 and rows[2, 5]["pin_sum"] == 5  # MKV KVL
    # 9 vs 1: both members, mean 5 k-mers, union of 6; 9 = member 1's sequence
    assert rows[9, 1]["shared"] == 5 and rows[9, 1]["n_members"] == 2
    assert rows[9, 1]["survival"] == pytest.approx(1.0)


def test_fit_alpha_recovers_the_shape() -> None:
    rng = np.random.default_rng(0)
    a = rng.uniform(0.6, 1.0, 5000)
    s = np.clip(survival(a, 11, 1.7) + rng.normal(0, 0.03, a.size), 0, 1)
    got, rmse = alpha.fit_alpha(a, s, 11)
    assert got == pytest.approx(1.7, rel=0.1) and rmse == pytest.approx(0.03, rel=0.2)
