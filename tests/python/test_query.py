"""Query counts against an index built from the fixture proteins."""

import itertools
import json
import time
from functools import partial
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from typer.testing import CliRunner

from kmer_functional_profiler import _core, query, reference
from kmer_functional_profiler.cli import app
from kmer_functional_profiler.index import (
    PIN_BITS,
    Index,
    IndexParams,
    PackedTable,
    UnitTable,
    build_index,
    write_unit_columns,
)
from kmer_functional_profiler.query import (
    UFIRST_SCORE,
    WTA_SCORE,
    Timer,
    _Summed,
    assign_best,
    em,
    em_pin,
    fit_present_prior,
    gather,
    posterior_zi,
    presence,
    profile,
    unit_hits,
)

DATA = Path(__file__).resolve().parents[1] / "data"
READS = (DATA / "reads_1.fastq.gz", DATA / "reads_2.fastq.gz")
K = 7


def proteins() -> dict[str, str]:
    records = (DATA / "proteins.faa").read_text().split(">")[1:]
    return {r.split("\n", 1)[0]: r.split("\n", 1)[1].replace("\n", "") for r in records}


@pytest.fixture(scope="module")
def members(tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("q") / "members.parquet"
    rows = [(i, i, True, seq) for i, seq in enumerate(proteins().values())]
    pl.DataFrame(
        rows, schema=["protein_id", "cluster_rep", "full_length", "sequence"], orient="row"
    ).write_parquet(path)
    return path


def build(members: Path, **params: object) -> Index:
    out = members.parent / "_".join(f"{k}{v}" for k, v in params.items())
    build_index(members, out, IndexParams(k=K, **params))  # type: ignore[arg-type]
    return Index.load(out)


def test_dense_counts_match_direct_intersection(members: Path) -> None:
    index = build(members, t_base=1.0, fp_bits=64)
    result = profile(index, *READS, all_estimators=True)
    hits = list(_core.FastxHits(*READS, k=K))
    hashes = [h for b in hits for h in b["hash"].tolist()]
    reads = [r for b in hits for r in b["read"].tolist()]
    expected = {}
    for unit, seq in enumerate(proteins().values()):
        kmers = set(reference.protein_kmers(seq.encode(), K))
        mine = [(h, r) for h, r in zip(hashes, reads, strict=True) if h in kmers]
        if mine:
            expected[unit] = (len(mine), len({h for h, _ in mine}), len({r for _, r in mine}))
    got = {row[0]: row[1:] for row in result.select("unit", "hits", "kmers_hit", "reads").rows()}
    assert got == expected
    # Gather hands every distinct hit k-mer to exactly one unit.
    all_kmers = set().union(*(reference.protein_kmers(s.encode(), K) for s in proteins().values()))
    assert result["kmers_unique"].sum() == len(set(hashes) & all_kmers)
    assert (result["kmers_unique"] <= result["kmers_hit"]).all()
    # EM keeps every hit: expected hits over the detected units equal the observed ones.
    expected_hits = sum(h in all_kmers for h in hashes)
    assert (result["coverage_em"] * result["m_g"]).sum() == pytest.approx(expected_hits)
    # So do the one-pass rules, which each give every hit k-mer to exactly one unit.
    for rule in ("wta", "ufirst"):
        assert result[f"kmers_{rule}"].sum() == result["kmers_unique"].sum()
        assert (result[f"coverage_{rule}"] * result["m_g"]).sum() == pytest.approx(expected_hits)
    # Every source protein is found by its reads.
    assert len(got) == len(proteins())


def test_sampled_index_hits_respect_unit_thresholds(members: Path) -> None:
    index = build(members, t_base=0.2, n_min=8)
    result = profile(index, *READS)
    exact = profile(build(members, t_base=0.2, n_min=8, fp_bits=64), *READS)
    # 16-bit fingerprints only add false hits (rare), never lose true ones.
    both = result.join(exact, on="unit", how="right", suffix="_exact")
    assert (both["hits"] >= both["hits_exact"]).all()
    assert (both["hits"] - both["hits_exact"]).sum() <= 2
    assert (result["kmers_hit"] <= result["m_g"]).all()
    assert result["containment"].is_between(0, 1).all()


def test_dense_tier_fits_em_on_every_kmer(members: Path) -> None:
    # A sparse tier 2 with a fully dense tier fits EM on the same k-mers as a dense index.
    both = build(members, t_base=0.2, n_min=0, t_dense=1.0, fp_bits=64)
    dense = build(members, t_base=1.0, fp_bits=64)
    assert both.dense is not None and dense.dense is None
    assert (both.units["m_dense"] == dense.units["m_g"]).all()
    assert both.units["m_g"].sum() < both.units["m_dense"].sum()
    got = profile(both, *READS, all_estimators=True)
    want = profile(dense, *READS, all_estimators=True)
    cols = ["unit", "coverage_em", "coverage_zi"]
    joined = got.filter(pl.col("kmers_unique") > 0).select(*cols, "kmers_dense")
    joined = joined.join(want.select(cols), on="unit", suffix="_want")
    assert joined.height == len(proteins())
    for col in cols[1:]:
        assert joined[col].to_list() == pytest.approx(joined[f"{col}_want"].to_list())
    assert (joined["kmers_dense"] > 0).all()


def test_copies_count_member_equivalents(tmp_path: Path) -> None:
    # Units of two unrelated proteins, both sequenced (a KO seen in two genomes): with every
    # kept k-mer present, copies = kept k-mers / an average member's = 2.
    seqs = list(proteins().values())
    rows = [(i, i // 2, True, seq) for i, seq in enumerate(seqs)]
    path = tmp_path / "pairs.parquet"
    schema = ["protein_id", "cluster_rep", "full_length", "sequence"]
    pl.DataFrame(rows, schema=schema, orient="row").write_parquet(path)
    build_index(path, tmp_path / "idx", IndexParams(k=K, t_base=1.0, fp_bits=64))
    result = profile(Index.load(tmp_path / "idx"), *READS, all_estimators=True)
    result = result.filter(pl.col("present_zi") == 1)
    assert result.height > 0
    assert result["copies_zi"].to_list() == pytest.approx([2.0] * result.height)
    assert (result["abundance_zi"] == 2 * result["coverage_zi"]).all()


def test_posterior_intervals_bracket_estimates(members: Path) -> None:
    index = build(members, t_base=1.0, fp_bits=64)
    plain = profile(index, *READS, all_estimators=True)
    got = profile(index, *READS, draws=60, all_estimators=True)
    assert "coverage_zi_lo" not in plain.columns
    # intervals and groups change nothing else
    assert got.drop("^.*_(lo|hi)$", "ambiguity_group", "group_size", "own_evidence").equals(plain)
    assert got["own_evidence"].drop_nulls().is_between(0, 1).all()
    found = got.filter(pl.col("coverage_zi") > 0)
    assert (found["coverage_zi_lo"] <= found["coverage_zi_hi"]).all()
    inside = found["coverage_zi"].is_between(found["coverage_zi_lo"], found["coverage_zi_hi"])
    assert inside.mean() >= 0.9  # type: ignore[operator]
    assert got.equals(profile(index, *READS, draws=60, all_estimators=True))  # seeded


def test_shared_evidence_groups_lopsided_pair() -> None:
    # Unit 0: 20 k-mers of its own and 3 shared with unit 1; unit 2 stands alone. 4 hits
    # (one read each) per k-mer. Unit 1 lives mostly on the shared k-mers, which barely move
    # unit 0's coverage, so the pair is lopsided.
    def run(own_of_1: list[int]) -> tuple[dict[int, int | None], dict[int, float]]:
        kmer_units = [(h, [0]) for h in range(20)] + [(h, [0, 1]) for h in (20, 21, 22)]
        kmer_units += [(h, [1]) for h in own_of_1] + [(h, [2]) for h in range(30, 40)]
        rows = [(u, h, h * 10 + r, 1) for h, us in kmer_units for u in us for r in range(4)]
        rows = [x for x in rows if x[1] not in own_of_1 or x[2] % 10 == 0]  # 1 stray hit
        schema = {"unit": pl.UInt32, "hash": pl.UInt64, "read": pl.UInt64, "n": pl.UInt32}
        hit_reads = pl.DataFrame(rows, schema=schema, orient="row")
        got = posterior_zi(hit_reads, np.array([25, 5, 12]), np.array([25.0, 5.0, 12.0]), 200)
        return (
            dict(got.select("unit", "ambiguity_group").iter_rows()),
            dict(got.select("unit", "own_evidence").iter_rows()),
        )

    # One k-mer of its own with a single stray hit (a typical false positive after gather):
    # mostly shared evidence.
    group, own = run([25])
    assert group[0] == group[1] == 0 and group[2] is None
    assert own[2] == 1 and 0.8 < own[0] < 1 and own[1] < 0.5
    # None of its own: explained away in every draw, still grouped with unit 0.
    group, own = run([])
    assert group[0] == group[1] == 0 and group[2] is None


def test_gather_explains_away_shared_kmers() -> None:
    kmers = pl.DataFrame(
        {"unit": [0, 0, 0, 0, 1, 1, 2, 2, 3, 3], "hash": [1, 2, 3, 4, 3, 4, 4, 5, 6, 7]},
        schema={"unit": pl.UInt32, "hash": pl.UInt64},
    )
    t_g = np.array([0.01, 0.01, 0.01, 0.01])
    got = {u: (n, r) for u, n, r in gather(kmers, t_g).iter_rows()}
    # 1 is a subset of 0: explained away. 2 keeps only hash 5; 3 shares nothing.
    assert got == {0: (4, 0), 2: (1, 2), 3: (2, 1)}
    # Sampled 10x more sparsely, unit 1's 2 k-mers stand for more than unit 0's 4.
    got = {u: n for u, n, _ in gather(kmers, np.array([0.01, 0.001, 0.01, 0.01])).iter_rows()}
    assert got == {1: 2, 0: 2, 2: 1, 3: 2}


def test_assign_best_rules() -> None:
    # Same k-mers as the gather test: unit 0 holds 1-4, unit 1 holds 3-4, unit 2 holds 4-5.
    kmers = pl.DataFrame(
        {"unit": [0, 0, 0, 0, 1, 1, 2, 2, 3, 3], "hash": [1, 2, 3, 4, 3, 4, 4, 5, 6, 7]},
        schema={"unit": pl.UInt32, "hash": pl.UInt64},
    )

    def won(m_g: list[int], t_g: list[float], score: pl.Expr) -> dict[int, list[int]]:
        rated = kmers.with_columns(
            m_g=pl.col("unit").replace_strict(range(4), m_g),
            t_g=pl.col("unit").replace_strict(range(4), t_g),
        )
        rows = assign_best(rated, score).group_by("unit").agg(pl.col("hash").sort())
        return dict(rows.sort("unit").iter_rows())

    m_g, t_g = [8, 2, 4, 2], [0.01] * 4
    # Winner-take-all: unit 1 is fully contained (2/2), so it takes 3 and 4 from unit 0 (4/8).
    assert won(m_g, t_g, WTA_SCORE) == {0: [1, 2], 1: [3, 4], 2: [5], 3: [6, 7]}
    # Uniqueness-first: unit 0 has the most specific support (1 + 1 + 1/2 + 1/3) and goes first.
    assert won(m_g, t_g, UFIRST_SCORE) == {0: [1, 2, 3, 4], 2: [5], 3: [6, 7]}
    # Sampled 10x more sparsely, unit 1's support (1/2 + 1/3) / 0.001 outranks unit 0's.
    assert won(m_g, [0.01, 0.001, 0.01, 0.01], UFIRST_SCORE) == {
        0: [1, 2],
        1: [3, 4],
        2: [5],
        3: [6, 7],
    }


def test_em_splits_shared_kmer_hits() -> None:
    # Unit 0 holds k-mers 1, 2 and unit 1 holds 2, 3; hits 10, 15, 5. The Poisson MLE
    # solves 10/l0 + 15/(l0+l1) = 2 = 5/l1 + 15/(l0+l1): l0 = 10, l1 = 5.
    kmers = pl.DataFrame(
        {"unit": [0, 0, 1, 1], "hash": [1, 2, 2, 3], "hits": [10, 15, 15, 5]},
        schema={"unit": pl.UInt32, "hash": pl.UInt64, "hits": pl.UInt32},
    )
    got = em(kmers, np.array([2, 2]))
    assert got["coverage"].to_list() == pytest.approx([10.0, 5.0], rel=1e-4)
    # Unhit kept k-mers pull coverage down: unit 1 with 4 kept k-mers.
    got = em(kmers.filter(pl.col("unit") == 1), np.array([2, 4]))
    assert got["coverage"].to_list() == pytest.approx([20 / 4])


def test_zero_inflated_em_fits_present_kmers() -> None:
    # 10 of 40 kept k-mers present, 3 hits each: the other 30 are structural zeros, not
    # low coverage. Zero-truncated Poisson: c / (1 - exp(-c)) = 3.
    kmers = pl.DataFrame(
        {"unit": [0] * 10, "hash": range(10), "hits": [3] * 10},
        schema={"unit": pl.UInt32, "hash": pl.UInt64, "hits": pl.UInt32},
    )
    got = em(kmers, np.array([40]), zero_inflated=True).row(0, named=True)
    c = got["coverage"]
    assert c / -np.expm1(-c) == pytest.approx(3.0, rel=1e-4)
    assert got["present"] == pytest.approx(10 / (40 * -np.expm1(-c)), rel=1e-4)
    assert em(kmers, np.array([40]))["coverage"].item() == pytest.approx(30 / 40)
    # All kept k-mers hit once: no excess zeros, so present stays 1 and it matches plain EM.
    got = em(kmers.with_columns(hits=pl.lit(1, pl.UInt32)), np.array([10]), zero_inflated=True)
    assert got.row(0) == (0, pytest.approx(1.0), 1.0)


def test_pin_presence_matches_zero_inflated_em_at_one_level() -> None:
    # Every kept k-mer at the same p_in level: presence cannot vary between k-mers, so the
    # fit is zero-inflated EM's (10 of 40 k-mers present, 3 hits each).
    kmers = pl.DataFrame(
        {"unit": [0] * 10, "hash": range(10), "hits": [3] * 10, "pin_q": [7] * 10},
        schema={"unit": pl.UInt32, "hash": pl.UInt64, "hits": pl.UInt32, "pin_q": pl.UInt8},
    )
    hist = np.zeros((1, 16))
    hist[0, 7] = 40
    got = em_pin(kmers, hist).row(0, named=True)
    want = em(kmers.drop("pin_q"), np.array([40]), zero_inflated=True).row(0, named=True)
    assert got["coverage"] == pytest.approx(want["coverage"], rel=1e-4)
    assert got["present"] == pytest.approx(want["present"], rel=1e-4)


def test_pin_presence_escapes_core_kmers() -> None:
    # Every kept k-mer core (p_in = 1) but only 10 of 40 present: presence must fall, as in
    # zero-inflated EM, rather than stay at 1 (it did when level 15 read as exactly 1).
    kmers = pl.DataFrame(
        {"unit": [0] * 10, "hash": range(10), "hits": [3] * 10, "pin_q": [15] * 10},
        schema={"unit": pl.UInt32, "hash": pl.UInt64, "hits": pl.UInt32, "pin_q": pl.UInt8},
    )
    hist = np.zeros((1, 16))
    hist[0, 15] = 40
    got = em_pin(kmers, hist).row(0, named=True)
    want = em(kmers.drop("pin_q"), np.array([40]), zero_inflated=True).row(0, named=True)
    assert got["coverage"] == pytest.approx(want["coverage"], rel=1e-3)
    assert got["present"] == pytest.approx(want["present"], rel=1e-3)


def test_pin_presence_splits_shared_kmer_by_p_in() -> None:
    # Units 0 and 1 have one mid-p_in k-mer each (5 hits) and share one (10 hits) that is
    # core (p_in ~1) in unit 0 and private (p_in ~0) in unit 1: unit 0 should take it.
    kmers = pl.DataFrame(
        {
            "unit": [0, 0, 1, 1],
            "hash": [1, 3, 2, 3],
            "hits": [5, 10, 5, 10],
            "pin_q": [8, 15, 8, 0],
        },
        schema={"unit": pl.UInt32, "hash": pl.UInt64, "hits": pl.UInt32, "pin_q": pl.UInt8},
    )
    hist = np.zeros((2, 16))
    hist[0, [8, 15]] = 1
    hist[1, [8, 0]] = 1
    got = dict(em_pin(kmers, hist).select("unit", "coverage").iter_rows())
    assert got[0] > 1.2 * got[1]  # unit 1's alpha rises to explain the k-mer, so not all
    flat = dict(
        em(kmers, np.array([2, 2]), zero_inflated=True).select("unit", "coverage").iter_rows()
    )
    assert flat[0] == pytest.approx(flat[1])


def test_present_prior_shrinks_thin_units() -> None:
    fit = pl.DataFrame({"coverage": [5.0] * 20, "present": [0.2, 0.4] * 10})
    a, b = fit_present_prior(fit)  # type: ignore[misc]
    assert a / (a + b) == pytest.approx(0.3)
    assert fit_present_prior(fit.head(5)) is None
    # One hit k-mer out of 100 at coverage ~1: present is barely measured, so the prior
    # pulls it towards 0.3 and coverage down to match; 1000 hit k-mers of 10,000 hold.
    for n, m, moved in ((1, 100, True), (1000, 10_000, False)):
        kmers = pl.DataFrame(
            {"unit": [0] * n, "hash": range(n), "hits": [2] * n},
            schema={"unit": pl.UInt32, "hash": pl.UInt64, "hits": pl.UInt32},
        )
        free = em(kmers, np.array([m]), zero_inflated=True).row(0, named=True)
        shrunk = em(kmers, np.array([m]), zero_inflated=True, prior=(a, b)).row(0, named=True)
        assert (abs(shrunk["present"] - free["present"]) > 0.05) == moved
        assert abs(shrunk["present"] - 0.3) <= abs(free["present"] - 0.3) + 1e-9


def reference_lookup(table: PackedTable, hashes: np.ndarray) -> np.ndarray:
    """The numpy lookup the Rust kernel replaced: materialised sorted keys, searchsorted."""
    keys = np.repeat(np.arange(len(table.offsets) - 1, dtype=np.uint64), np.diff(table.offsets))
    keys = keys << np.uint64(table.fp_bits) | table.fingerprints.astype(np.uint64)
    query = hashes >> np.uint64(table.shift)
    i = np.minimum(np.searchsorted(keys, query), len(keys) - 1)
    found = (keys[i] == query) & (hashes <= np.uint64(table.max_hash))
    return np.where(found, table.set_ids[i].astype(np.int64), -1)


def reference_unit_hits(
    table: PackedTable, max_hash_g: np.ndarray, hashes: np.ndarray, reads: np.ndarray
) -> pl.DataFrame:
    """The numpy set expansion and ``max_hash_g`` check the Rust kernel replaced."""
    set_ids = reference_lookup(table, hashes)
    found = set_ids >= 0
    hashes, reads, set_ids = hashes[found], reads[found], set_ids[found]
    starts = table.set_offsets[set_ids].astype(np.int64)
    lengths = table.set_offsets[set_ids + 1].astype(np.int64) - starts
    within = np.arange(lengths.sum()) - np.repeat(np.cumsum(lengths) - lengths, lengths)
    values = table.set_values[np.repeat(starts, lengths) + within].astype(np.uint64)
    units = (values >> np.uint64(PIN_BITS)).astype(np.uint32)
    pin_q = (values & np.uint64(2**PIN_BITS - 1)).astype(np.uint8)
    hit = np.repeat(np.arange(len(hashes)), lengths)
    keep = hashes[hit] <= max_hash_g[units]
    holders = np.bincount(hit[keep], minlength=len(hashes)).astype(np.uint32)
    hit = hit[keep]
    return pl.DataFrame(
        {
            "unit": units[keep],
            "hash": hashes[hit],
            "read": reads[hit],
            "pin_q": pin_q[keep],
            "holders": holders[hit],
        }
    )


def test_unit_hits_matches_numpy_reference(tmp_path: Path) -> None:
    # 4-bit fingerprints: many keys share a bucket and fingerprint, so sets are unions.
    rng = np.random.default_rng(5)
    max_hash, n_units = 2**58 - 1, 300
    hashes = rng.integers(0, max_hash, 20_000, dtype=np.uint64)
    units = rng.integers(0, n_units, len(hashes), dtype=np.uint64)
    values = units << np.uint64(PIN_BITS) | rng.integers(0, 16, len(hashes), dtype=np.uint64)
    table = PackedTable.build(hashes, values, max_hash, 4)
    max_hash_g = rng.integers(0, max_hash, n_units, dtype=np.uint64)
    query = np.concatenate([hashes[::3], rng.integers(0, 2**59, 5000, dtype=np.uint64)])
    reads = np.arange(len(query), dtype=np.uint64)
    want = reference_unit_hits(table, max_hash_g, query, reads)
    assert want.height > 1000 and (want["holders"] > 1).any()
    assert unit_hits(table, max_hash_g, query, reads).equals(want)
    # Memory-mapped, as a loaded index reads it.
    layout = table.save(tmp_path, "t")
    mapped = PackedTable.load(tmp_path, "t", layout)
    assert isinstance(mapped.set_values, np.memmap)
    assert unit_hits(mapped, max_hash_g, query, reads).equals(want)
    assert np.array_equal(table.lookup(query), reference_lookup(table, query))


def test_cli_query(members: Path, tmp_path: Path) -> None:
    runner = CliRunner()
    idx = tmp_path / "idx"
    assert runner.invoke(app, ["index", str(members), str(idx), "--k", str(K)]).exit_code == 0
    out = tmp_path / "p.tsv"
    stats = tmp_path / "stats.json"
    args = ["query", str(idx), *map(str, READS), "--out", str(out), "--stats", str(stats)]
    result = runner.invoke(app, [*args, "--draws", "3"])
    assert result.exit_code == 0, result.output
    table = pl.read_csv(out, separator="\t")
    assert {"cluster_rep", "hits", "containment", "coverage"} <= set(table.columns)
    got = json.loads(stats.read_text())
    assert {"load", "hash", "lookup", "gather", "fit_zi", "posterior", "total"} <= set(
        got["stages"]
    )
    counts = got["counts"]
    assert counts["hit_units"] == table.height
    # Units are renumbered for the model; the output carries the index's ids again.
    units = pl.read_parquet(idx / "units.parquet").select("unit", "cluster_rep")
    joined = table.join(units, on="unit", suffix="_index")
    assert joined.height == table.height
    assert (joined["cluster_rep"] == joined["cluster_rep_index"]).all()
    assert counts["hit_rows"] == table["hits"].sum()
    assert counts["sampled_kmers"] >= counts["hit_kmers"] > 0
    assert counts["largest_component_units"] <= counts["hit_units"]
    assert got["running"] == ""
    # Tiers read into memory give the same profile as memory-mapped ones.
    mapped, loaded = tmp_path / "mapped.tsv", tmp_path / "loaded.tsv"
    base = ["query", str(idx), *map(str, READS), "--out"]
    assert runner.invoke(app, [*base, str(mapped)]).exit_code == 0
    assert runner.invoke(app, [*base, str(loaded), "--in-memory"]).exit_code == 0
    assert pl.read_csv(loaded, separator="\t").equals(pl.read_csv(mapped, separator="\t"))


def test_timer_writes_stats_mid_stage(tmp_path: Path) -> None:
    # A query killed out of memory leaves the stages done and the one it died in.
    stats = tmp_path / "stats.json"
    timer = Timer(stats)
    with timer("load"), timer("keys"):
        got = json.loads(stats.read_text())
    assert got["running"] == "keys" and got["stages"] == {}
    assert json.loads(stats.read_text())["stages"].keys() == {"load", "keys"}


def test_presence_doubts_few_and_shared_hits() -> None:
    # 400 units with 5-30 own k-mers, 100 with one unique own k-mer, 20 with one k-mer 40
    # index units hold, 50 with two. 2000 index units at t_g 0.01 and 3 M reads expect
    # ~48 background units with one k-mer.
    rng = np.random.default_rng(0)
    h = np.array([*rng.integers(5, 30, 400), *[1] * 120, *[2] * 50])
    holders = np.ones(h.sum(), dtype=np.uint32)
    holders[np.cumsum(h)[400 + 100 : 400 + 120] - 1] = 40
    own = pl.DataFrame(
        {"unit": np.repeat(np.arange(len(h)), h), "hash": np.arange(h.sum()), "holders": holders},
        schema={"unit": pl.UInt32, "hash": pl.UInt64, "holders": pl.UInt32},
    )
    t_g = np.full(len(h), 0.01)
    got = presence(own, t_g, 3_000_000, 2000)["present_prob"].to_numpy()
    many, unique, shared, two = got[:400], got[400], got[500], got[520]
    assert many.min() > 0.95
    assert shared < unique < two < 1
    # Ten times fewer reads: less background, so one k-mer is more credible.
    fewer = presence(own, t_g, 300_000, 2000)["present_prob"].to_numpy()
    assert fewer[400] > unique


def test_posterior_draws_ignore_input_order() -> None:
    # Seeded draws must not depend on row order, which Polars group_by does not fix.
    rows = [(u, h, h * 10 + r, 1 + r % 2) for u in (0, 1) for h in range(u * 3, u * 3 + 8)
            for r in range(3)]  # fmt: skip
    schema = {"unit": pl.UInt32, "hash": pl.UInt64, "read": pl.UInt64, "n": pl.UInt32}
    hit_reads = pl.DataFrame(rows, schema=schema, orient="row")
    m_g, pin_sum = np.array([10, 10]), np.array([10.0, 10.0])
    got = posterior_zi(hit_reads, m_g, pin_sum, 20)
    assert got.equals(
        posterior_zi(hit_reads.sample(fraction=1.0, shuffle=True, seed=1), m_g, pin_sum, 20)
    )


def test_unit_columns_match_parquet(members: Path, tmp_path: Path) -> None:
    # Memory-mapped unit columns give the parquet's rows and the same profile as an index
    # written before them (columns read from units.parquet).
    index = build(members, t_base=0.2, n_min=0, t_dense=1.0)
    old = tmp_path / "old"
    old.mkdir()
    for f in Path(index.units.path).parent.iterdir():
        if not f.name.startswith("units."):
            (old / f.name).write_bytes(f.read_bytes())
    (old / "units.parquet").write_bytes(index.units.path.read_bytes())
    assert isinstance(index.units["m_g"], np.memmap)
    some = np.array([0, 2, 3], dtype=np.uint32)
    assert index.units.rows(some).equals(index.units.frame()[some.tolist()])
    got = profile(index, *READS, draws=20)
    want = profile(Index.load(old), *READS, draws=20)
    assert got.equals(want)


def test_unit_columns_read_in_slices(
    members: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Columns too long for one Polars array are read in slices: same arrays either way.
    index = build(members, t_base=0.2, n_min=0)
    want = {n: np.asarray(index.units[n]) for n in ("pin_hist", "t_g", "m_g")}
    monkeypatch.setattr("kmer_functional_profiler.index.COLUMN_SLICE", 3)
    (tmp_path / "units.parquet").write_bytes(index.units.path.read_bytes())
    fallback = UnitTable(tmp_path)
    for n, a in want.items():
        assert np.array_equal(fallback[n], a)
    write_unit_columns(tmp_path)
    for n, a in want.items():
        assert np.array_equal(np.load(tmp_path / f"units.{n}.npy"), a)


def test_em_fits_components_independently() -> None:
    # Components share no k-mer, so fitting them together gives each one's own fit: a
    # slow component (few hits, many kept k-mers) no longer holds a fast one to its pace.
    rng = np.random.default_rng(3)
    n, pairs = 60, 400
    kmers = pl.DataFrame(
        {
            "unit": rng.integers(0, n, pairs),
            "hash": rng.integers(0, 300, pairs),
            "hits": rng.integers(1, 6, pairs),
            "pin_q": rng.integers(0, 16, pairs),
        },
        schema={"unit": pl.UInt32, "hash": pl.UInt64, "hits": pl.UInt32, "pin_q": pl.UInt8},
    ).unique(["unit", "hash"], keep="first", maintain_order=True)
    kmers = kmers.with_columns(pl.col("hits").max().over("hash"))
    hist = rng.integers(0, 4, (n, 16)).astype(np.float64) + 1
    shifted = kmers.with_columns(unit=pl.col("unit") + n, hash=pl.col("hash") + 1000)
    together = pl.concat([kmers, shifted.head(50)])
    m_g = hist.sum(axis=1)
    for fit, a, b in (
        (em_pin, (together, np.vstack([hist, hist])), (kmers, hist)),
        (em, (together, np.r_[m_g, m_g]), (kmers, m_g)),
        (partial(em, zero_inflated=True), (together, np.r_[m_g, m_g]), (kmers, m_g)),
    ):
        got, alone = fit(*a).filter(pl.col("unit") < n), fit(*b)
        assert got["unit"].equals(alone["unit"])
        for col in ("coverage", "present"):
            assert got[col].to_list() == pytest.approx(alone[col].to_list(), rel=1e-12)


def test_default_fits_em_only(members: Path) -> None:
    # The benchmark estimators are opt-in and leave the shipped columns unchanged.
    index = build(members, t_base=1.0, fp_bits=64)
    default, full = profile(index, *READS), profile(index, *READS, all_estimators=True)
    assert not any(c.endswith(("_zi", "_zib", "_zip", "_wta", "_ufirst")) for c in default.columns)
    assert {"coverage_zi", "coverage_zib", "coverage_zip", "kmers_wta"} <= set(full.columns)
    assert default.equals(full.select(default.columns))


def test_summed_batches_equal_one_aggregation() -> None:
    rng = np.random.default_rng(0)
    frame = pl.DataFrame({"unit": rng.integers(0, 5, 300), "n": rng.integers(1, 4, 300)})
    summed = _Summed(["unit"], {"n": pl.col("n").sum()})
    summed.min_rows = 1  # merge at every batch
    for start in range(0, 300, 40):
        summed.add(frame[start : start + 40])
    got, want = summed.total().sort("unit"), frame.group_by("unit").agg(pl.col("n").sum())
    assert got.equals(want.sort("unit"))


def test_timer_keeps_the_sampled_anonymous_peak(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = itertools.count()
    # A spike seen only by the sampler mid-stage (call 0 is Timer's own start).
    monkeypatch.setattr(query, "anon_rss", lambda: 100 if next(calls) == 2 else 1)
    timer = Timer()
    with timer("stage"):
        time.sleep(0.3)
    assert timer.stages["stage"]["peak_anon"] == 100
    monkeypatch.setattr(query, "anon_rss", lambda: None)  # no /proc (macOS)
    timer = Timer()
    with timer("stage"):
        pass
    assert "peak_anon" not in timer.stages["stage"]
