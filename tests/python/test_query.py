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
    sample_summary,
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


@pytest.fixture(scope="module")
def shared(members: Path) -> Path:
    """The fixture proteins plus a copy of each, every 20th residue changed, as a cluster of
    its own: units share k-mers, so components have two units."""
    path = members.parent / "shared" / "members.parquet"
    path.parent.mkdir(exist_ok=True)
    seqs = list(proteins().values())
    copies = ["".join("W" if i % 20 == 10 else c for i, c in enumerate(q)) for q in seqs]
    rows = [(i, i, True, q) for i, q in enumerate(seqs + copies)]
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


def test_copies_finite_when_kept_kmers_are_fragment_only(tmp_path: Path) -> None:
    # The unit's full-length member is too short for a k-mer, so every kept k-mer comes from
    # its fragment and has p_in 0: pin_sum was 0 and copies_zi infinite.
    rows = [(0, 0, True, "MK"), (1, 0, False, max(proteins().values(), key=len))]
    path = tmp_path / "frag.parquet"
    schema = ["protein_id", "cluster_rep", "full_length", "sequence"]
    pl.DataFrame(rows, schema=schema, orient="row").write_parquet(path)
    build_index(path, tmp_path / "idx", IndexParams(k=K, t_base=1.0, fp_bits=64))
    got = profile(Index.load(tmp_path / "idx"), *READS, all_estimators=True, with_aai=True)
    assert got.height == 1
    assert got.select(pl.col("copies_zi", "abundance_zi").is_finite().all()).row(0) == (True, True)
    assert got["aai"].drop_nulls().is_between(0, 1).all()


def test_posterior_intervals_bracket_estimates(members: Path) -> None:
    index = build(members, t_base=1.0, fp_bits=64)
    plain = profile(index, *READS, all_estimators=True)
    got = profile(index, *READS, draws=60, all_estimators=True)
    assert "coverage_zi_lo" not in plain.columns
    # intervals and groups change nothing else
    posterior = ("^coverage_zi_(lo|hi)$", "^.*abundance_zi_(lo|hi)$", "^group_coverage.*$")
    assert got.drop(*posterior, "ambiguity_group", "group_size", "own_evidence").equals(plain)
    assert got["own_evidence"].drop_nulls().is_between(0, 1).all()
    found = got.filter(pl.col("coverage_zi") > 0)
    assert (found["coverage_zi_lo"] <= found["coverage_zi_hi"]).all()
    inside = found["coverage_zi"].is_between(found["coverage_zi_lo"], found["coverage_zi_hi"])
    assert inside.mean() >= 0.9  # type: ignore[operator]
    assert got.equals(profile(index, *READS, draws=60, all_estimators=True))  # seeded


def test_ztp_lambda_inverts_the_truncated_mean() -> None:
    lam = np.array([0.05, 0.5, 2.0, 10.0, 40.0])
    mean = lam / -np.expm1(-lam)
    assert query.ztp_lambda(mean) == pytest.approx(lam, rel=1e-9)
    assert (query.ztp_lambda(np.array([0.0, 1.0])) == 0).all()


def test_aai_naive_equals_aai_without_shared_kmers(members: Path) -> None:
    index = build(members, t_base=1.0, fp_bits=64)
    got = profile(index, *READS, with_aai=True)
    assert {"aai", "aai_naive", "aai_naive_lower_bound", "component", "hits_em"} <= set(got.columns)
    assert got["component"].n_unique() == got.height  # unrelated proteins: no links
    both = got.filter(pl.col("aai").is_not_null() & ~pl.col("aai_naive_lower_bound"))
    assert both.height > 0
    assert both["aai"].to_list() == pytest.approx(both["aai_naive"].to_list(), abs=1e-6)
    assert got["aai_naive"].is_between(0, 1).all()
    # Reads come from the indexed proteins themselves: identity ~1 where coverage allows.
    assert both["aai"].median() > 0.95  # type: ignore[operator]
    # min_aai drops rows; the plain profile has no aai but has aai_naive.
    assert profile(index, *READS, min_aai=1.1).height == 0
    assert "aai" not in profile(index, *READS).columns


def test_aai_of_a_unit_explained_away(tmp_path: Path) -> None:
    # Unit 1 is the first half of unit 0's protein: all its k-mers are unit 0's too.
    first = max(proteins().values(), key=len)
    rows = [(0, 0, True, first), (1, 1, True, first[: len(first) // 2])]
    path = tmp_path / "m.parquet"
    schema = ["protein_id", "cluster_rep", "full_length", "sequence"]
    pl.DataFrame(rows, schema=schema, orient="row").write_parquet(path)
    build_index(path, tmp_path / "idx", IndexParams(k=K, t_base=1.0, fp_bits=64))
    got = profile(Index.load(tmp_path / "idx"), *READS, with_aai=True).sort("unit")
    assert got.height == 2
    assert got["component"].to_list() == [0, 0]
    half = got.row(1, named=True)
    assert half["kmers_unique"] == 0  # gather gives its k-mers to unit 0
    assert half["aai"] is None
    assert half["aai_naive"] > 0.7  # its own hits still show it close (~1x: uncorrected)


def test_aai_interval_without_draws(members: Path) -> None:
    index = build(members, t_base=1.0, fp_bits=64)
    got = profile(index, *READS, with_aai=True).filter(pl.col("aai").is_not_null())
    assert got.height > 0
    assert (got["aai_lo"] <= got["aai"] + 1e-12).all() and (
        got["aai"] <= got["aai_hi"] + 1e-12
    ).all()
    assert got["aai_hi"].max() <= 1.0  # type: ignore[operator]


def test_aai_interval_survival_variance() -> None:
    # Substitutions at identity a over L positions: the share of the L - k + 1 windows that
    # survive has variance c (1 - c) d / L, d the overlap factor aai_interval uses.
    rng = np.random.default_rng(0)
    k, length, a = 11, 400, 0.93
    hit = rng.random((4000, length)) < a
    windows = np.lib.stride_tricks.sliding_window_view(hit, k, axis=1).all(axis=2)
    c = windows.mean(axis=1)
    n = length - k + 1
    lo, hi = query.aai_interval(
        coverage=np.array([1e6]), present=np.array([a**k]), m=np.array([n]),
        pin_sum=np.array([n]), n_kmers=np.array([n]), k=k,
    )  # fmt: skip
    # The interval of a from the survival term alone holds ~95% of simulated strains.
    inside = ((c ** (1 / k) >= lo[0]) & (c ** (1 / k) <= hi[0])).mean()
    assert 0.92 <= inside <= 0.98
    # More kept windows, narrower; sparse sampling (t << 1) removes the overlap term.
    args = {"coverage": np.full(3, 1e6), "present": np.full(3, a**k), "k": k}
    lo, hi = query.aai_interval(m=np.array([20, 200, 200]), pin_sum=np.array([20, 200, 200]),
                                n_kmers=np.array([20, 200, 2000]), **args)  # fmt: skip
    width = hi - lo
    assert width[0] > width[1] > width[2]


def write_members(path: Path, seqs: list[str], first_rep: int = 0) -> Path:
    rows = [(first_rep + i, first_rep + i, True, q) for i, q in enumerate(seqs)]
    schema = ["protein_id", "cluster_rep", "full_length", "sequence"]
    pl.DataFrame(rows, schema=schema, orient="row").write_parquet(path)
    return path


@pytest.mark.parametrize(("t_base", "t_dense"), [(1.0, 0.0), (0.3, 1.0)])
def test_joint_query_equals_one_index_of_the_union(
    tmp_path: Path, t_base: float, t_dense: float
) -> None:
    seqs = list(proteins().values())
    half = len(seqs) // 2
    params = IndexParams(k=K, t_base=t_base, t_dense=t_dense, fp_bits=64)
    for name, part, first in (("a", seqs[:half], 0), ("b", seqs[half:], half), ("ab", seqs, 0)):
        build_index(
            write_members(tmp_path / f"{name}.parquet", part, first), tmp_path / name, params
        )
    a, b, ab = (Index.load(tmp_path / n) for n in ("a", "b", "ab"))
    union = profile(ab, *READS, draws=5, with_aai=True)
    joint = profile(a, *READS, draws=5, with_aai=True, extra=[b])
    assert joint["source"].to_list() == [int(u >= a.units.height) for u in joint["unit"]]
    assert joint.drop("source").equals(union.select(joint.drop("source").columns))
    # An unrelated extra index leaves the first index's rows as they were, but for
    # present_prob, whose prior counts every index's units.
    unrelated = tmp_path / "unrelated"
    rng = np.random.default_rng(0)
    random_seqs = ["".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 300)) for _ in range(5)]
    build_index(write_members(tmp_path / "r.parquet", random_seqs), unrelated, params)
    alone = profile(a, *READS)
    with_unrelated = profile(a, *READS, extra=[Index.load(unrelated)])
    assert with_unrelated.filter(pl.col("source") == 1).height == 0
    keep = [c for c in alone.columns if c != "present_prob"]
    assert with_unrelated.select(keep).equals(alone.select(keep))


def test_joint_query_rejects_a_different_scheme(tmp_path: Path) -> None:
    seqs = list(proteins().values())[:3]
    path = write_members(tmp_path / "m.parquet", seqs)
    build_index(path, tmp_path / "k7", IndexParams(k=7, t_base=1.0))
    build_index(path, tmp_path / "k8", IndexParams(k=8, t_base=1.0))
    with pytest.raises(ValueError, match="differ in k"):
        profile(Index.load(tmp_path / "k7"), *READS, extra=[Index.load(tmp_path / "k8")])


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
        got = posterior_zi(
            *by_hash(hit_reads), np.array([25, 5, 12]), np.array([25.0, 5.0, 12.0]), 200
        )
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
    summary = tmp_path / "summary.json"
    result = runner.invoke(app, [*args, "--draws", "3", "--aai", "--summary", str(summary)])
    assert result.exit_code == 0, result.output
    table = pl.read_csv(out, separator="\t")
    assert {"cluster_rep", "hits", "containment", "coverage", "aai_lo", "component"} <= set(
        table.columns
    )
    sample = json.loads(summary.read_text())
    assert sample["explained_fraction"] > 0 and 0 <= sample["census_containment"] <= 1
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


def test_frame_modes_and_quality_mask_reach_the_profile(members: Path, tmp_path: Path) -> None:
    runner = CliRunner()
    idx = tmp_path / "idx"
    build_args = ["index", str(members), str(idx), "--k", str(K), "--t-base", "1"]
    assert runner.invoke(app, build_args).exit_code == 0
    hits = {}
    for name, extra in {
        "stopfree": [],
        "edges": ["--frames", "edges:5"],
        "masked": ["--min-qual", "41"],  # above every fixture quality: everything is N
        "joint": ["--extra-index", str(idx)],  # itself: every unit twice, each half the hits
    }.items():
        out = tmp_path / f"{name}.tsv"
        args = ["query", str(idx), *map(str, READS), "--out", str(out), *extra]
        assert runner.invoke(app, args).exit_code == 0
        hits[name] = pl.read_csv(out, separator="\t").select("unit", name=pl.col("hits"))
    # Edges hash a superset of the stop-free k-mers, so no unit loses hits.
    both = hits["stopfree"].join(hits["edges"], on="unit", how="left", suffix="_e")
    assert (both["name_e"] >= both["name"]).all()
    assert hits["edges"]["name"].sum() > hits["stopfree"]["name"].sum()
    assert hits["masked"].height == 0
    assert hits["joint"].height == 2 * hits["stopfree"].height


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
    report: dict[str, int] = {}
    presence(own, t_g, 3_000_000, 2000, report=report)
    assert report["presence_converged"] == 1 and 1 < report["presence_iterations"] < 500
    presence(own, t_g, 3_000_000, 2000, max_iter=1, report=report)
    assert report == {"presence_iterations": 1, "presence_converged": 0}


def by_hash(hit_reads: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """``posterior_zi``'s inputs from per-(unit, hash, read) rows."""
    return hit_reads.select("unit", "hash").unique(), hit_reads.select("hash", "read", "n").unique()


def test_posterior_draws_ignore_input_order() -> None:
    # Seeded draws must not depend on row order, which Polars group_by does not fix.
    rows = [(u, h, h * 10 + r, 1 + r % 2) for u in (0, 1) for h in range(u * 3, u * 3 + 8)
            for r in range(3)]  # fmt: skip
    schema = {"unit": pl.UInt32, "hash": pl.UInt64, "read": pl.UInt64, "n": pl.UInt32}
    hit_reads = pl.DataFrame(rows, schema=schema, orient="row")
    m_g, pin_sum = np.array([10, 10]), np.array([10.0, 10.0])
    got = posterior_zi(*by_hash(hit_reads), m_g, pin_sum, 20)
    assert got.equals(
        posterior_zi(
            *by_hash(hit_reads.sample(fraction=1.0, shuffle=True, seed=1)), m_g, pin_sum, 20
        )
    )


def test_low_memory_rereads_to_the_same_profile(shared: Path) -> None:
    index = build(shared, t_base=0.2, n_min=0)
    timer = Timer()
    got = profile(index, *READS, draws=5, batch_reads=7, low_memory=True, timer=timer)
    assert got.equals(profile(index, *READS, draws=5, batch_reads=7))
    assert "reread" in timer.stages and timer.counts["read_rows"] == 0


def test_read_weights_are_poisson_and_keyed_by_read() -> None:
    reads = np.arange(200_000, dtype=np.uint64)
    w = query._poisson1(3, 7, reads)
    assert abs(w.mean() - 1) < 0.01 and abs(w.var() - 1) < 0.02
    assert (query._poisson1(3, 7, reads[::-1]) == w[::-1]).all()  # a read's own weight
    assert not (query._poisson1(3, 8, reads) == w).all()  # a new draw, new weights


def test_posterior_component_alone_or_batched(monkeypatch: pytest.MonkeyPatch) -> None:
    # One component per batch: units 0 and 1 (sharing k-mers) get the same intervals whether
    # or not unit 2, hit by the same reads on k-mers of its own, is queried with them.
    monkeypatch.setattr(query, "POSTERIOR_BATCH_BYTES", 1)
    kmer_units = [(h, [0]) for h in range(6)] + [(h, [0, 1]) for h in (6, 7)]
    kmer_units += [(h, [2]) for h in range(10, 14)]
    rows = [(u, h, r, 1 + r % 2) for h, us in kmer_units for u in us for r in range(h % 3, 9, 2)]
    schema = {"unit": pl.UInt32, "hash": pl.UInt64, "read": pl.UInt64, "n": pl.UInt32}
    hit_reads = pl.DataFrame(rows, schema=schema, orient="row")
    m_g, pin_sum = np.array([10, 5, 8]), np.array([10.0, 5.0, 8.0])
    both = posterior_zi(*by_hash(hit_reads), m_g, pin_sum, 30)
    alone = posterior_zi(*by_hash(hit_reads.filter(pl.col("unit") < 2)), m_g, pin_sum, 30)
    assert both["group_size"].max() == 2
    assert both.filter(pl.col("unit") < 2).equals(alone)


def test_posterior_batches_in_parallel_as_in_series(monkeypatch: pytest.MonkeyPatch) -> None:
    # Batches are seeded by their own units, so worker processes give the serial result.
    monkeypatch.setattr(query, "POSTERIOR_BATCH_BYTES", 1)
    kmer_units = [(h, [u]) for u in range(5) for h in range(u * 10, u * 10 + 6)]
    kmer_units += [(h, [0, 1]) for h in (50, 51)]
    rows = [(u, h, r, 1 + r % 2) for h, us in kmer_units for u in us for r in range(h % 3, 9, 2)]
    schema = {"unit": pl.UInt32, "hash": pl.UInt64, "read": pl.UInt64, "n": pl.UInt32}
    hit_reads = pl.DataFrame(rows, schema=schema, orient="row")
    m_g, pin_sum = np.full(5, 10), np.full(5, 10.0)
    serial = posterior_zi(*by_hash(hit_reads), m_g, pin_sum, 30, workers=1)
    assert posterior_zi(*by_hash(hit_reads), m_g, pin_sum, 30, workers=3).equals(serial)


def test_em_reports_unconverged_components() -> None:
    # Two units sharing every k-mer split them slowly; with 3 iterations both are cut off.
    kmers = pl.DataFrame(
        {"unit": [0, 0, 1, 1, 2], "hash": [1, 2, 1, 2, 3], "hits": [5, 7, 5, 7, 4]},
        schema={"unit": pl.UInt32, "hash": pl.UInt64, "hits": pl.UInt32},
    )
    report: dict[str, int] = {}
    em(kmers, np.array([4, 6, 2]), max_iter=3, report=report)
    assert report == {"em_iterations": 3, "em_unconverged_units": 2}
    em(kmers, np.array([4, 6, 2]), report=report)
    assert report["em_unconverged_units"] == 2 and report["em_iterations"] > 3


def test_component_batches_change_nothing(shared: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # One component per batch: every fit gives what it gives on all units at once.
    indexes = [build(shared, t_base=0.2, n_min=0), build(shared, t_base=0.2, n_min=0, t_dense=1.0)]
    timer = Timer()
    profile(indexes[0], *READS, timer=timer)
    assert timer.counts["components"] < timer.counts["hit_units"]
    assert timer.counts["fit_batches"] == 1 and timer.counts["em_unconverged_units"] == 0
    assert 0 < timer.counts["em_iterations"] <= 1000
    assert timer.counts["presence_converged"] == 1
    assert timer.counts["links"] >= timer.counts["cut_links_0.2"] >= timer.counts["cut_links_0.005"]
    whole = [profile(i, *READS, all_estimators=True, draws=3) for i in indexes]
    monkeypatch.setattr(query, "MAX_BATCH_PAIRS", 1)
    for index, expected in zip(indexes, whole, strict=True):
        assert profile(index, *READS, all_estimators=True, draws=3).equals(expected)
    pairs = pl.DataFrame(
        {"unit": [0, 1, 1, 2, 3, 4], "hash": [10, 10, 11, 12, 11, 13]},
        schema={"unit": pl.UInt32, "hash": pl.UInt64},
    )
    batches = [sorted(b["unit"].unique().to_list()) for b in query.component_batches(pairs)]
    assert sorted(batches) == [[0, 1, 3], [2], [4]]


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
    assert not list(tmp_path.glob("*.partial"))


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
    monkeypatch.setattr(query, "anon_rss", lambda: None)  # neither /proc nor macOS
    timer = Timer()
    with timer("stage"):
        pass
    assert "peak_anon" not in timer.stages["stage"]


def _families(n: int, seed: int = 0, chained: bool = False) -> tuple[pl.DataFrame, np.ndarray]:
    """Families of units sharing core k-mers (slow for EM to split), Poisson hits; with
    ``chained``, consecutive units also share a k-mer, so all form one component."""
    rng = np.random.default_rng(seed)
    unit, kmer, m_g, u, h = [], [], [], 0, 0
    for _ in range(n):
        core = int(rng.integers(5, 40))
        for j in range(int(rng.integers(1, 6))):
            own = int(rng.integers(0, 20))
            ks = [*range(h, h + core), *range(h + 100 * (j + 1), h + 100 * (j + 1) + own)]
            unit += [u] * len(ks)
            kmer += ks
            m_g.append(len(ks) + int(rng.integers(0, 10)))
            u += 1
        h += 1000
    if chained:
        unit += [*range(u - 1), *range(1, u)]
        kmer += [10**9 + i for i in range(u - 1)] * 2
    pairs = pl.DataFrame(
        {"unit": unit, "hash": kmer}, schema={"unit": pl.UInt32, "hash": pl.UInt64}
    )
    cov = rng.gamma(0.5, 2.0, u) * (rng.random(u) < 0.6)
    per_hash = (
        pairs.with_columns(c=pl.Series(cov[unit]))
        .group_by("hash")
        .agg(pl.col("c").sum())
        .sort("hash")
    )
    per_hash = per_hash.with_columns(
        hits=pl.Series(rng.poisson(per_hash["c"].to_numpy()) + chained)
    )
    kmers = pairs.join(per_hash.select("hash", hits=pl.col("hits").cast(pl.UInt32)), on="hash")
    return kmers.filter(pl.col("hits") > 0), np.array(m_g)


def _plain_two_steps(step, la, se, *_):  # type: ignore[no-untyped-def]
    l1, s1, _att = step(la, se)
    l2, s2, att = step(l1, s1)
    return l1, s1, l2, s2, att


def test_squarem_reaches_plain_em_fixed_point(monkeypatch: pytest.MonkeyPatch) -> None:
    # SQUAREM: the same coverages as plain EM run to convergence, in fewer EM steps.
    kmers, m_g = _families(200)
    fast: dict[str, int] = {}
    got = em(kmers, m_g, report=fast, max_iter=100_000)
    monkeypatch.setattr(query, "_squarem", _plain_two_steps)
    slow: dict[str, int] = {}
    want = em(kmers, m_g, report=slow, max_iter=100_000)
    assert fast["em_unconverged_units"] == slow["em_unconverged_units"] == 0
    assert fast["em_iterations"] * 2 < slow["em_iterations"]
    assert np.allclose(got["coverage"], want["coverage"], rtol=0, atol=1e-5)


def test_blockwise_em_matches_whole(monkeypatch: pytest.MonkeyPatch) -> None:
    # One component fitted in blocks of ~300 pairs reaches the whole fit's coverages (plain EM
    # has one fixed point). The zero-inflated fits may reach another of their fixed points.
    kmers, m_g = _families(60, chained=True)
    assert query.components(kmers)[0] == 1
    whole = em(kmers, m_g, tol=1e-10)
    monkeypatch.setattr(query, "MAX_FIT_PAIRS", 300)
    report: dict[str, int] = {}
    blocks = em(kmers, m_g, tol=1e-10, report=report)
    assert report["em_block_rounds"] > 1 and report["em_unconverged_units"] == 0
    assert np.allclose(blocks["coverage"], whole["coverage"], rtol=0, atol=1e-6)
    zi = em(kmers, m_g, zero_inflated=True)
    assert zi.height == whole.height and np.isfinite(zi["coverage"].to_numpy()).all()


def test_link_cuts_on_largest_component() -> None:
    # Units 0-1 share 4 of 4 k-mers (strength 1); 1-2 share one k-mer of 10 (0.1, 1 hit);
    # unit 3 is a component of its own.
    rows = [(0, h, 3) for h in range(4)] + [(1, h, 3) for h in [*range(4), *range(20, 30)]]
    rows += [(1, 4, 1), (2, 4, 1)] + [(2, h, 2) for h in range(5, 14)] + [(3, 99, 5)]
    kmers = pl.DataFrame(
        rows, schema={"unit": pl.UInt32, "hash": pl.UInt64, "hits": pl.UInt32}, orient="row"
    )
    got = query.link_cuts(kmers)
    assert got["links"] == 2 and got["links_one_kmer"] == 1 and got["links_le2_hits"] == 1
    assert got["detected_largest_component_units"] == 3
    assert got["detected_largest_component_pairs"] == 4 + 15 + 10
    assert got["cut_links_0.1"] == 0 and got["cut_largest_units_0.1"] == 3
    assert got["cut_links_0.2"] == 1 and got["cut_largest_units_0.2"] == 2


def test_posterior_in_blocks(monkeypatch: pytest.MonkeyPatch) -> None:
    # Reads reweighted in chunks: identical (integer sums). Sweeps in k-mer blocks: the same
    # model, other random numbers; the same units, groups and similar intervals.
    kmer_units = [(h, [0]) for h in range(6)] + [(h, [0, 1]) for h in range(6, 12)]
    kmer_units += [(h, [1, 2]) for h in range(12, 16)] + [(h, [2]) for h in range(16, 20)]
    rows = [(u, h, r, 1 + r % 2) for h, us in kmer_units for u in us for r in range(h % 3, 9, 2)]
    schema = {"unit": pl.UInt32, "hash": pl.UInt64, "read": pl.UInt64, "n": pl.UInt32}
    hit_reads = pl.DataFrame(rows, schema=schema, orient="row")
    m_g, pin_sum = np.array([10, 14, 8]), np.array([10.0, 14.0, 8.0])
    want = posterior_zi(*by_hash(hit_reads), m_g, pin_sum, 200)
    monkeypatch.setattr(query, "READ_CHUNK", 5)
    assert posterior_zi(*by_hash(hit_reads), m_g, pin_sum, 200).equals(want)
    monkeypatch.setattr(query, "SWEEP_BLOCK_PAIRS", 4)
    got = posterior_zi(*by_hash(hit_reads), m_g, pin_sum, 200)
    assert got.select("unit", "ambiguity_group", "group_size").equals(
        want.select("unit", "ambiguity_group", "group_size")
    )
    assert np.allclose(got["coverage_zi_hi"], want["coverage_zi_hi"], rtol=0.3)


CODONS = [
    a + b + c
    for a in "ACGT"
    for b in "ACGT"
    for c in "ACGT"
    if a + b + c not in ("TAA", "TAG", "TGA")
]


def test_sample_summary_by_hand() -> None:
    result = pl.DataFrame({"coverage_zi": [2.0, 0.0], "len_mean": [100.0, 50.0]})
    stats = {"bases": 10_000.0, "lost": 0.0, "expected_errors": 0.0}
    got = sample_summary(result, 11, 100, 1, stats, (200, 150), 0.2)
    # 100 reads of 100 bp: 98 / 3 whole codons per frame, 98 / 3 - 10 windows; no errors.
    rho = (100 / 3) / (98 / 3 - 10)
    assert got["error_thinning"] == 1.0
    assert got["explained_bases"] == pytest.approx(2.0 * 300 * rho)
    assert got["explained_fraction"] == pytest.approx(600 * rho / 10_000)
    assert got["unknown_fraction"] == pytest.approx(1 - 600 * rho / 10_000)
    assert got["census_containment"] == pytest.approx(0.75)
    assert got["census_frame_miss"] == pytest.approx(0.8 ** (98 / 3 - 10))
    # Errors thin k-mer depth: per-base loss d gives r = (1 - d)^(3k), which both estimates
    # divide out; masked bases count in full, Phred errors at SENSE_CHANGE.
    noisy = sample_summary(
        result, 11, 100, 1, stats | {"lost": 20.0, "expected_errors": 40.0}, (200, 150), 0.2
    )
    r = (1 - (20 + query.SENSE_CHANGE * 40) / 10_000) ** 33
    assert noisy["error_thinning"] == pytest.approx(r)
    assert noisy["explained_fraction"] == pytest.approx(got["explained_fraction"] / r)  # type: ignore[operator]
    assert noisy["census_containment"] == pytest.approx(min(1.0, 0.75 / r))
    # Without base counts (sourmash hashing) or a census, those values are None.
    bare = sample_summary(result, 11, 100, 1, None, None, 0.2)
    assert bare["explained_fraction"] is None and bare["census_containment"] is None


def test_unknown_fraction_of_known_and_random_reads(tmp_path: Path) -> None:
    # Five stop-free genes of 400 codons, each indexed as its own unit; half of the reads
    # come from them (either strand), half are random DNA.
    rng = np.random.default_rng(11)
    genes = ["".join(rng.choice(CODONS, 400)) for _ in range(5)]
    table = str.maketrans("ACGT", "TGCA")
    members = tmp_path / "members.parquet"
    pl.DataFrame(
        [(i, i, True, _core.translate_frames(g.encode())[0].decode()) for i, g in enumerate(genes)],
        schema=["protein_id", "cluster_rep", "full_length", "sequence"],
        orient="row",
    ).write_parquet(members)
    build_index(members, tmp_path / "idx", IndexParams(k=K, t_base=0.05))
    index = Index.load(tmp_path / "idx")

    def reads(n: int, source: list[str] | None) -> list[str]:
        out = []
        for _ in range(n):
            if source is None:
                out.append("".join(rng.choice(list("ACGT"), 150)))
            else:
                g = source[rng.integers(len(source))]
                start = rng.integers(len(g) - 150 + 1)
                read = g[start : start + 150]
                out.append(read if rng.random() < 0.5 else read.translate(table)[::-1])
        return out

    def summary(seqs: list[str]) -> dict[str, float | int | None]:
        path = tmp_path / "reads.fa"
        path.write_text("".join(f">{i}\n{q}\n" for i, q in enumerate(seqs)))
        got: dict[str, float | int | None] = {}
        profile(index, path, summary=got)
        return got

    known, unknown = reads(400, genes), reads(400, None)
    mixed = summary(known + unknown)
    # Reads never cross the genes' ends here, so end k-mers are thinly covered, and taken as
    # absent they lift coverage_zi: explained runs ~10% high on known reads.
    assert mixed["bases"] == 800 * 150 and mixed["error_thinning"] == 1.0
    assert mixed["explained_fraction"] == pytest.approx(0.5, abs=0.08)
    only_known = summary(known)
    assert only_known["explained_fraction"] == pytest.approx(1.0, abs=0.12)
    assert only_known["census_containment"] > 0.85  # type: ignore[operator]
    only_random = summary(unknown)
    assert only_random["census_containment"] < 0.05  # type: ignore[operator]
    assert (only_random["explained_fraction"] or 0.0) < 0.05
