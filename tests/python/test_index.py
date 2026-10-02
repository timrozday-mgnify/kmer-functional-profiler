"""Index build on hand-built clusters, against the analytical sample size, and lookup."""

import json
import random
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from typer.testing import CliRunner

from kmer_functional_profiler import index, reference
from kmer_functional_profiler.cli import app
from kmer_functional_profiler.index import PIN_BITS, Index, IndexParams, PackedTable, build_index

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
K = 6
RNG = random.Random(7)
X, D, Y, P, Q, S = ("".join(RNG.choices(AMINO_ACIDS, k=n)) for n in (30, 30, 30, 20, 20, 40))
# (protein_id, cluster_rep, full_length, sequence): unit A = 1, B = 10, singleton S = 20.
HAND = [
    (1, 1, True, X + D),
    (2, 1, True, X + D),
    (3, 1, True, X + D),
    (4, 1, True, X + D + P),
    (5, 1, False, Q + X[:10]),
    (10, 10, True, D + Y),
    (11, 10, True, D + Y),
    (20, 20, True, S),
]


def write_members(path: Path, rows: list[tuple[int, int, bool, str]]) -> Path:
    pl.DataFrame(
        rows, schema=["protein_id", "cluster_rep", "full_length", "sequence"], orient="row"
    ).write_parquet(path)
    return path


def kmers(seq: str) -> set[int]:
    return set(reference.protein_kmers(seq.encode(), K))


@pytest.fixture
def hand(tmp_path: Path) -> tuple[Index, pl.DataFrame]:
    members = write_members(tmp_path / "members.parquet", HAND)
    build_index(
        members, tmp_path / "idx", IndexParams(k=K, t_base=1.0, n_min=1), postings_parquet=True
    )
    index = Index.load(tmp_path / "idx")
    postings = pl.read_parquet(tmp_path / "idx" / "postings.parquet").join(
        index.units.frame().select("unit", "cluster_rep"), on="unit"
    )
    return index, postings


def rows_for(postings: pl.DataFrame, rep: int, part: str) -> pl.DataFrame:
    return postings.filter(pl.col("cluster_rep") == rep, pl.col("hash").is_in(kmers(part)))


def test_scores_on_hand_built_clusters(hand: tuple[Index, pl.DataFrame]) -> None:
    _, postings = hand
    core, private, partial_only, shared = (rows_for(postings, 1, part) for part in (X, P, Q, D))
    assert core.height == len(kmers(X))
    assert (core["p_in"] == 1).all() and (core["score"] == 0).all()
    assert (private["p_in"] == 0.25).all() and (private["score"] == -2).all()
    # half a member of 4, not 0
    assert (partial_only["p_in"] == 0.125).all() and (partial_only["score"] == -3).all()
    assert (shared["n_groups"] == 2).all() and (shared["score"] == -1).all()
    assert (rows_for(postings, 10, D)["score"] == -1).all()


def test_unit_table(hand: tuple[Index, pl.DataFrame]) -> None:
    index, _ = hand
    units = {row["cluster_rep"]: row for row in index.units.frame().iter_rows(named=True)}
    a, b, s = units[1], units[10], units[20]
    assert a["n_members"] == 5
    assert a["n_kmers"] == len(set().union(*(kmers(r[3]) for r in HAND if r[1] == 1)))
    # Kept k-mers per counting member: 3 x (X + D), 1 x (X + D + P); their mean is pin_sum,
    # plus half a member's share (0.5 / 4) of each k-mer only the partial member 5 holds.
    held = np.array([len(kmers(X + D))] * 3 + [len(kmers(X + D + P))])
    partial_only = kmers(Q + X[:10]) - kmers(X + D + P)
    assert a["pin_sum"] == pytest.approx(held.mean() + 0.125 * len(partial_only))
    assert a["len_cv"] == pytest.approx(held.std() / held.mean())
    assert b["len_cv"] == s["len_cv"] == 0


def test_promiscuous_kmers_are_dropped(tmp_path: Path) -> None:
    # Unit 30 holds only D, shared with A and B, so it has no posting and is not indexed.
    members = write_members(tmp_path / "members.parquet", [*HAND, (30, 30, True, D)])
    params = IndexParams(k=K, t_base=1.0, n_min=1, max_groups=2)
    stats = build_index(members, tmp_path / "idx", params)
    assert stats["promiscuous_dropped"] == 3 * len(kmers(D))
    assert (stats["n_clusters"], stats["n_units"]) == (4, 3)
    units = Index.load(tmp_path / "idx").units.frame()
    assert units["cluster_rep"].to_list() == [1, 10, 20]
    assert units["unit"].to_list() == [0, 1, 2]
    a = (kmers(X + D + P) | kmers(Q + X[:10])) - kmers(D)
    assert units["m_g"].to_list() == [len(a), len(kmers(D + Y) - kmers(D)), len(kmers(S))]


def test_tier2_lookup_returns_postings(hand: tuple[Index, pl.DataFrame]) -> None:
    index, postings = hand
    set_ids = index.tier2.lookup(postings["hash"].to_numpy())
    assert (set_ids >= 0).all()
    for set_id, unit, pin_q in zip(set_ids, postings["unit"], postings["pin_q"], strict=True):
        assert (unit << PIN_BITS | pin_q) in index.tier2.values(int(set_id))


def test_candidates_match_analytical_expectation(tmp_path: Path) -> None:
    rng = random.Random(1)
    rows = []
    for rep in range(1, 400):
        base = "".join(rng.choices(AMINO_ACIDS, k=rng.randint(30, 400)))
        for m in range(rng.choice([1, 1, 1, 2, 5])):
            mutated = "".join(rng.choice(AMINO_ACIDS) if rng.random() < 0.05 else a for a in base)
            rows.append((rep * 100 + m, rep * 100, rng.random() < 0.7, mutated))
    members = write_members(tmp_path / "members.parquet", rows)
    params = IndexParams(k=K, t_base=0.02, n_min=8)
    stats = build_index(members, tmp_path / "a", params, postings_parquet=True)
    expected = stats["candidates_expected"]
    assert isinstance(expected, float)
    assert abs(stats["candidates"] - expected) < 5 * expected**0.5  # type: ignore[operator]

    units = Index.load(tmp_path / "a").units.frame()
    raised = params.oversample * params.n_min / units["n_kmers"]
    expected_t = raised.clip(upper_bound=params.t_cap).clip(lower_bound=params.t_base)
    t_g = pl.Series(np.where(units["n_members"] > 1, expected_t, params.t_base))
    assert (units["t_g"] - t_g).abs().max() < 1e-12  # type: ignore[operator]
    floored = units.filter(pl.col("t_g") > params.t_base)
    assert floored.height > 0 and (floored["m_g"] <= params.n_min).all()
    postings = pl.read_parquet(tmp_path / "a" / "postings.parquet").join(units, on="unit")
    assert (postings["hash"] <= postings["max_hash_g"]).all()

    # Unit-aligned batching does not change the result.
    build_index(
        members,
        tmp_path / "b",
        IndexParams(k=K, t_base=0.02, n_min=8, batch_residues=900),
        postings_parquet=True,
    )
    for name in ("postings", "units"):
        a, b = (pl.read_parquet(tmp_path / d / f"{name}.parquet") for d in "ab")
        assert a.equals(b)


def test_packed_table_false_hits_and_range() -> None:
    rng = np.random.default_rng(3)
    max_hash = 2**60 - 1
    hashes = rng.integers(0, max_hash, 5000, dtype=np.uint64)
    table = PackedTable.build(hashes, np.arange(5000, dtype=np.uint64), max_hash, 16)
    assert table.lead_bits == 4
    assert (table.lookup(hashes) >= 0).all()
    other = rng.integers(0, max_hash, 100_000, dtype=np.uint64)
    assert (table.lookup(other) >= 0).mean() < 1e-3
    assert table.lookup(np.array([max_hash + 1], dtype=np.uint64))[0] == -1


def test_cli_index(tmp_path: Path) -> None:
    members = write_members(tmp_path / "members.parquet", HAND)
    result = CliRunner().invoke(
        app, ["index", str(members), str(tmp_path / "idx"), "--k", str(K), "--t-base", "1"]
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["n_units"] == 3


def test_floored_units_keep_best_scoring_kmers(tmp_path: Path) -> None:
    members = write_members(tmp_path / "members.parquet", HAND)
    # Every k-mer is a candidate (t_cap = 1); floored units keep their n_min best.
    params = IndexParams(k=K, t_base=0.001, n_min=5, t_cap=1.0, oversample=1000)
    build_index(members, tmp_path / "idx", params, postings_parquet=True)
    index = Index.load(tmp_path / "idx")
    units = index.units.frame().sort("cluster_rep")
    assert units["m_g"].to_list()[:2] == [5, 5]
    postings = pl.read_parquet(tmp_path / "idx" / "postings.parquet").join(
        units.select("unit", "cluster_rep"), on="unit"
    )
    assert set(postings.filter(pl.col("cluster_rep") == 1)["hash"]) <= kmers(X + D) - kmers(D)
    assert (postings["score"] == 0).all()


def test_adapter_peptides_are_masked(tmp_path: Path) -> None:
    # Translated Nextera read-through (frame 0 of CTGTCTCTTATACACATCTCCGAGCCCACGAGAC).
    adapter = "LSLIHISEPTRPLY"
    rows = [(1, 1, True, X + adapter), (2, 2, True, Y + adapter), (3, 3, True, S)]
    members = write_members(tmp_path / "members.parquet", rows)
    for mask in (True, False):
        # At k = 11 no k-mer keeps 6 unmasked adapter residues (the mask width).
        params = IndexParams(k=11, t_base=1.0, n_min=1, mask_adapters=mask)
        stats = build_index(members, tmp_path / str(mask), params, postings_parquet=True)
        hashes = set(pl.read_parquet(tmp_path / str(mask) / "postings.parquet")["hash"])
        adapter_kmers = set(reference.protein_kmers(adapter.encode(), 11))
        assert bool(hashes & adapter_kmers) is not mask
        assert stats["n_adapter_masked_proteins"] == (2 if mask else 0)


def test_packed_table_layout_is_compact() -> None:
    rng = np.random.default_rng(5)
    hashes = rng.integers(0, 2**64 - 1, 50_000, dtype=np.uint64)
    table = PackedTable.build(hashes, np.arange(50_000, dtype=np.uint64) % 1000, 2**64 - 1, 16)
    assert 2 <= len(hashes) / 2**table.bucket_bits <= 4
    assert table.set_values.dtype == np.uint16
    assert table.nbytes() / len(hashes) < 10


def test_madvise_random_reaches_memory_mapped_arrays(tmp_path: Path) -> None:
    # The hook acts only when a loaded array's base is the mmap; a numpy change there
    # would skip it silently.
    np.save(tmp_path / "a.npy", np.arange(10))
    mapped = np.load(tmp_path / "a.npy", mmap_mode="r")
    assert isinstance(mapped.base, index.mmap.mmap)
    assert index._madvise_random(mapped) is mapped
