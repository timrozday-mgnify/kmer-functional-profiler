"""The aai-model workflow's survival counting and fit."""

import json
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


def test_pin_sum_falls_back_to_all_members_when_p_was_the_only_full_length_one() -> None:
    # Cluster 1: P = 2 (full-length) and fragment 3. Without P no full-length member is left:
    # pin_sum is the fragment's k-mers, not 0 / 0 (step 25's run lost such pairs as NaN).
    members = pl.DataFrame(
        {
            "protein_id": [1, 2, 3],
            "cluster_rep": [1, 1, 1],
            "full_length": [False, True, False],
            "sequence": ["MKVLA", "MKVLAWE", "MKVLAW"],
        }
    )
    pairs = pl.DataFrame({"protein_id": [2], "cluster_rep": [1], "identity": [0.9]})
    got = fit_model.pair_survival(pairs, members, K, "protein").row(0, named=True)
    assert got["n_members"] == 2 and got["pin_sum"] == pytest.approx((3 + 4) / 2)
    assert got["shared"] == 4 and got["n_windows"] == 5  # MKV KVL VLA LAW of 5 windows


def test_fit_model_recovers_markov_beta_and_end_loss() -> None:
    # Pairs drawn from the model's own process, the last 1 - c of each P's windows outside
    # the aligned block: the binned fit finds phi, region and c.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from test_query import strains, windows_of

    rng = np.random.default_rng(4)
    true, c = SurvivalModel(11, None, 8.0, categories=8, concentration=5.0), 0.93
    rows = []
    for t in np.linspace(0.02, 0.5, 40):
        same = strains(true, float(t), 60, 300, rng)
        win = windows_of(same, 11)
        win[:, int(c * win.shape[1]) :] = False
        for i in range(len(win)):
            w = win[i]
            both = {f"both_{j}": (w[:-j] & w[j:]).mean() for j in fit_model.LAGS}
            rows.append({"identity": same[i].mean(), "window_survival": w.mean(),
                         "n_windows": len(w)} | both)  # fmt: skip
    got = fit_model.fit_model(pl.DataFrame(rows), 11, categories=8)
    assert got["ends"] == pytest.approx(c, abs=0.02)
    assert got["concentration"] == pytest.approx(5.0, rel=0.4)
    assert got["region"] == pytest.approx(8.0, rel=0.4)
    assert got["rmse"] < 0.3 * got["rmse_independent"]


def test_pairs_take_the_nearest_member_of_each_candidate_cluster(tmp_path: Path) -> None:
    # P = 7 against cluster 1 (members 1, 2, 7) and cluster 5 (member 5). Nearest in 1 is
    # member 2 (P itself and the short-coverage hit to 1 at 0.99 do not count).
    cands = tmp_path / "candidates.parquet"
    pl.DataFrame(
        {"protein_id": [7, 7], "cluster_rep": [1, 5], "identity_rep": [0.70, 0.65]}
    ).write_parquet(cands)
    clusters = tmp_path / "member_clusters.parquet"
    pl.DataFrame({"protein_id": [1, 2, 7, 5], "cluster_rep": [1, 1, 1, 5]}).write_parquet(clusters)
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        "7\t7\t100.0\t100\t100\t100\t200\n"  # P itself
        "7\t1\t99.0\t50\t100\t100\t90\n"  # coverage 0.5
        "7\t1\t70.0\t100\t100\t100\t120\n"
        "7\t2\t85.0\t95\t100\t100\t150\n"
        "7\t5\t65.0\t100\t100\t100\t100\n"
    )
    out = tmp_path / "pairs.parquet"
    args = ["pairs", "--hits", str(hits), "--candidates", str(cands),
            "--member-clusters", str(clusters), "--out", str(out),
            "--stats-out", str(tmp_path / "stats.json")]  # fmt: skip
    sys.argv = ["fit_model.py", *args]
    fit_model.main()
    got = {r["cluster_rep"]: r for r in pl.read_parquet(out).iter_rows(named=True)}
    assert got[1]["nearest"] == 2 and got[1]["identity"] == pytest.approx(0.85)
    assert got[1]["identity_rep"] == pytest.approx(0.70)
    assert got[5]["nearest"] == 5 and got[5]["identity"] == pytest.approx(0.65)
    stats = json.loads((tmp_path / "stats.json").read_text())
    assert stats["candidates"] == 2 and stats["with_member_hit"] == 2 and stats["pairs"] == 2
    assert stats["with_member_hit_by_size"] == {"1": 1.0, "2-3": 1.0}


def test_union_fit_recovers_the_union_term_and_fit_writes_both_stages(tmp_path: Path) -> None:
    # Multi-member pairs whose window survival is the union model's plus noise: union_fit
    # finds (g0, g1). fit() then writes the one-member model and the union term.
    rng = np.random.default_rng(7)
    base = SurvivalModel(11, None, 5.0, categories=8, concentration=4.0, ends=0.95)
    true = SurvivalModel(11, None, 5.0, categories=8, concentration=4.0, ends=0.95,
                         union=0.08, union_slope=0.3)  # fmt: skip
    a = rng.uniform(0.6, 1.0, 20_000)
    n = np.exp(rng.uniform(np.log(2), np.log(500), a.size)).round()
    noisy = true.union_survival(a, n) + rng.normal(0, 0.05, a.size)
    many = pl.DataFrame({"identity": a, "n_members": n, "window_survival": noisy})
    got = fit_model.union_fit(many, base, 11)
    assert got["union"] == pytest.approx(0.08, abs=0.02)
    assert got["union_slope"] == pytest.approx(0.3, abs=0.1)
    assert got["rmse"] < 0.01 < got["rmse_independent"]
    # end to end: one-member pairs from the base model's process, plus the union pairs
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from test_query import strains, windows_of

    rows = []
    for t in np.linspace(0.02, 0.5, 30):
        same = strains(base, float(t), 40, 300, rng)
        win = windows_of(same, 11)
        win[:, int(0.95 * win.shape[1]) :] = False
        for i in range(len(win)):
            w = win[i]
            both = {f"both_{j}": (w[:-j] & w[j:]).mean() for j in fit_model.LAGS}
            rows.append({"identity": same[i].mean(), "window_survival": w.mean(),
                         "n_windows": len(w), "n_members": 1.0} | both)  # fmt: skip
    one = pl.DataFrame(rows)
    path = tmp_path / "survival.parquet"
    pl.concat([one, many], how="diagonal").write_parquet(path)
    out, strata = tmp_path / "aai_model.json", tmp_path / "strata.tsv"
    sys.argv = ["fit_model.py", "fit", "--survival", str(path), "--out", str(out),
                "--strata-out", str(strata), "--min-stratum", "500"]  # fmt: skip
    fit_model.main()
    model = __import__("json").loads(out.read_text())
    assert model["survival"] == "markov_beta" and model["ends"] == pytest.approx(0.95, abs=0.02)
    assert model["union"] > 0 and model["union_slope"] > 0
    names = pl.read_csv(strata, separator="\t")["stratum"].to_list()
    assert names[0] == "members 1" and "members 2+" in names and "members 101+" in names
