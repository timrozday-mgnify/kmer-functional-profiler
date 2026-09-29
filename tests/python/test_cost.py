"""Cost model: statistics pass against the build, and predicted against built sizes."""

import random
from pathlib import Path

import polars as pl

from kmer_functional_profiler.cost import cluster_stats, predict_cost
from kmer_functional_profiler.index import Index, IndexParams, build_index

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def test_prediction_matches_build(tmp_path: Path) -> None:
    rng = random.Random(5)
    rows = []
    for rep in range(1, 800):
        base = "".join(rng.choices(AMINO_ACIDS, k=rng.randint(20, 600)))
        for m in range(rng.choice([1, 1, 1, 2, 3, 8])):
            seq = "".join(rng.choice(AMINO_ACIDS) if rng.random() < 0.05 else a for a in base)
            rows.append((rep * 100 + m, rep * 100, rng.random() < 0.7, seq))
    members = tmp_path / "members.parquet"
    pl.DataFrame(
        rows, schema=["protein_id", "cluster_rep", "full_length", "sequence"], orient="row"
    ).write_parquet(members)
    params = IndexParams(k=6, t_base=0.01, n_min=8, t_dense=0.1)
    stats = build_index(members, tmp_path / "idx", params)
    units = Index.load(tmp_path / "idx").units

    clusters, pairs = cluster_stats(members, params, rate=0.05)
    assert clusters.select("cluster_rep", "n_members", "n_kmers").equals(
        units.select("cluster_rep", "n_members", "n_kmers")
    )
    assert pairs.height > 0 and pairs["hash"].max() <= 0.05 * 2**64  # type: ignore[operator]

    # Random k-mers are rarely shared, so keep = 1 and one key per posting hold; value sets
    # repeat within a unit, which the build measures.
    sets = stats["tier2_sets"] / stats["postings"]  # type: ignore[operator]
    cost = predict_cost(clusters.lazy(), params, sets=sets)
    assert cost["units"] == stats["n_units"] and cost["floored"] == stats["n_floored"]
    assert cost["t_max"] == stats["t_max"]
    for predicted, actual in (
        (cost["postings"], stats["postings"]),
        (cost["dense"], stats["dense_postings"]),
        (cost["tier2_bytes"], stats["tier2_bytes"]),
    ):
        assert abs(predicted / actual - 1) < 0.05, (predicted, actual)  # type: ignore[operator]
    # Dense sets repeat at a different rate, so its bytes are only bounded.
    assert cost["dense_bytes"] >= stats["dense_bytes"]  # type: ignore[operator]
