"""Query counts against an index built from the fixture proteins."""

from pathlib import Path

import numpy as np
import polars as pl
import pytest
from typer.testing import CliRunner

from kmer_functional_profiler import _core, reference
from kmer_functional_profiler.cli import app
from kmer_functional_profiler.index import Index, IndexParams, build_index
from kmer_functional_profiler.query import gather, profile

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


def test_cli_query(members: Path, tmp_path: Path) -> None:
    runner = CliRunner()
    idx = tmp_path / "idx"
    assert runner.invoke(app, ["index", str(members), str(idx), "--k", str(K)]).exit_code == 0
    out = tmp_path / "p.tsv"
    result = runner.invoke(app, ["query", str(idx), *map(str, READS), "--out", str(out)])
    assert result.exit_code == 0, result.output
    table = pl.read_csv(out, separator="\t")
    assert {"cluster_rep", "hits", "containment", "coverage"} <= set(table.columns)
