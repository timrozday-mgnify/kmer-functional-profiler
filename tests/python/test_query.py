"""Query counts against an index built from the fixture proteins."""

from pathlib import Path

import numpy as np
import polars as pl
import pytest
from typer.testing import CliRunner

from kmer_functional_profiler import _core, reference
from kmer_functional_profiler.cli import app
from kmer_functional_profiler.index import Index, IndexParams, build_index
from kmer_functional_profiler.query import (
    UFIRST_SCORE,
    WTA_SCORE,
    assign_best,
    em,
    em_pin,
    fit_present_prior,
    gather,
    posterior_zi,
    presence,
    profile,
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
    result = profile(index, *READS)
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
    got, want = profile(both, *READS), profile(dense, *READS)
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
    result = profile(Index.load(tmp_path / "idx"), *READS).filter(pl.col("present_zi") == 1)
    assert result.height > 0
    assert result["copies_zi"].to_list() == pytest.approx([2.0] * result.height)
    assert (result["abundance_zi"] == 2 * result["coverage_zi"]).all()


def test_posterior_intervals_bracket_estimates(members: Path) -> None:
    index = build(members, t_base=1.0, fp_bits=64)
    plain = profile(index, *READS)
    got = profile(index, *READS, draws=60)
    assert "coverage_zi_lo" not in plain.columns
    # intervals and groups change nothing else
    assert got.drop("^.*_(lo|hi)$", "ambiguity_group", "group_size", "own_evidence").equals(plain)
    assert got["own_evidence"].drop_nulls().is_between(0, 1).all()
    found = got.filter(pl.col("coverage_zi") > 0)
    assert (found["coverage_zi_lo"] <= found["coverage_zi_hi"]).all()
    inside = found["coverage_zi"].is_between(found["coverage_zi_lo"], found["coverage_zi_hi"])
    assert inside.mean() >= 0.9  # type: ignore[operator]
    assert got.equals(profile(index, *READS, draws=60))  # seeded


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


def test_cli_query(members: Path, tmp_path: Path) -> None:
    runner = CliRunner()
    idx = tmp_path / "idx"
    assert runner.invoke(app, ["index", str(members), str(idx), "--k", str(K)]).exit_code == 0
    out = tmp_path / "p.tsv"
    result = runner.invoke(app, ["query", str(idx), *map(str, READS), "--out", str(out)])
    assert result.exit_code == 0, result.output
    table = pl.read_csv(out, separator="\t")
    assert {"cluster_rep", "hits", "containment", "coverage"} <= set(table.columns)


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
