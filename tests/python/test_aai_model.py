"""The aai-model workflow's survival counting and fit."""

import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from kmer_functional_profiler.survival import SurvivalModel

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflows" / "aai-model"))
import fit_model  # noqa: E402

K = 3


def test_pair_survival_leaves_the_query_out_of_its_own_cluster() -> None:
    # Cluster 1: members 1 (MKVLAWE), 2 (MKVLAWQ); cluster 5: member 5 (MKVLRST). P = 2 in
    # its own cluster 1 and against cluster 5; P = 9 (member 1's sequence) against 1.
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
    got = fit_model.pair_survival(pairs, members, K, "protein")
    rows = {(r["protein_id"], r["cluster_rep"]): r for r in got.iter_rows(named=True)}
    # 2 vs 1 without itself: member 1 alone, 5 k-mers; 2 shares MKV KVL VLA LAW (4)
    assert rows[2, 1]["shared"] == 4 and rows[2, 1]["pin_sum"] == 5 and rows[2, 1]["n_members"] == 1
    assert rows[2, 5]["shared"] == 2 and rows[2, 5]["pin_sum"] == 5  # MKV KVL
    # 9 vs 1: both members, mean 5 k-mers, union of 6; 9 = member 1's sequence
    assert rows[9, 1]["shared"] == 5 and rows[9, 1]["n_members"] == 2
    assert rows[9, 1]["survival"] == pytest.approx(1.0)
    # 2's windows MKV KVL VLA LAW AWQ against member 1: found 1 1 1 1 0
    assert rows[2, 1]["window_survival"] == pytest.approx(0.8)
    assert rows[2, 1]["both_1"] == pytest.approx(3 / 4)  # (0,1) (1,2) (2,3) of 4 pairs
    assert rows[2, 1]["both_3"] == pytest.approx(1 / 2)  # (0,3) of (0,3) (1,4)
    assert rows[2, 1]["both_6"] is None  # no windows 6 apart


def test_fit_model_recovers_shape_and_region() -> None:
    # Pairs drawn from the model's own generative process: the fit finds its parameters.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from test_query import strains, windows_of

    rng = np.random.default_rng(4)
    true = SurvivalModel(11, 1.0, 30.0, categories=8)
    rows = []
    for t in np.linspace(0.02, 0.5, 25):
        same = strains(true, float(t), 40, 400, rng)
        win = windows_of(same, 11)
        for i in range(len(win)):
            w = win[i]
            both = {f"both_{j}": (w[:-j] & w[j:]).mean() for j in fit_model.LAGS}
            rows.append({"identity": same[i].mean(), "survival": w.mean()} | both)
    got = fit_model.fit_model(pl.DataFrame(rows), 11, categories=8)
    assert got["shape"] == pytest.approx(1.0, rel=0.3)
    assert got["region"] == pytest.approx(30.0, rel=0.3)
    assert got["rmse"] < 0.6 * got["rmse_independent"]  # per-pair noise dominates the rest
