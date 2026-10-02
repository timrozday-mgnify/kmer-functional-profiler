"""Mask sidecar: an index built from proteins a fixture "genome" encodes is masked by it."""

from pathlib import Path

import numpy as np
import polars as pl
import pytest
from typer.testing import CliRunner

from kmer_functional_profiler.cli import app
from kmer_functional_profiler.index import Index, IndexParams, build_index
from kmer_functional_profiler.mask import Mask, build_mask, windows
from kmer_functional_profiler.query import profile

DATA = Path(__file__).resolve().parents[1] / "data"
READS = (DATA / "reads_1.fastq.gz", DATA / "reads_2.fastq.gz")
K = 7


def records(path: Path) -> list[str]:
    return [r.split("\n", 1)[1].replace("\n", "") for r in path.read_text().split(">")[1:]]


def test_windows_cover_every_stretch(tmp_path: Path) -> None:
    rng = np.random.default_rng(0)
    seqs = ["".join(rng.choice(list("ACGT"), n)) for n in (95, 3, 40)]
    fasta = tmp_path / "g.fa"
    fasta.write_text(
        "".join(
            f">c{i}\n" + "\n".join(s[j : j + 7] for j in range(0, len(s), 7)) + "\n"
            for i, s in enumerate(seqs)
        )
    )
    got = [w.decode() for w in windows(fasta, 20, 5)]
    assert all(len(w) <= 20 for w in got)
    for s in seqs:
        for i in range(len(s) - 5):
            assert any(s[i : i + 6] in w for w in got)
    assert not any(seqs[0][-3:] + seqs[1] in w for w in got)  # records are never joined


@pytest.fixture(scope="module", params=[0.0, 1.0], ids=["tier2", "dense"])
def masked(
    request: pytest.FixtureRequest, tmp_path_factory: pytest.TempPathFactory
) -> tuple[Path, Path]:
    tmp = tmp_path_factory.mktemp("mask")
    seqs, cds = records(DATA / "proteins.faa"), records(DATA / "cds.fna")
    rows = [(i, i, True, q) for i, q in enumerate(seqs)]
    schema = ["protein_id", "cluster_rep", "full_length", "sequence"]
    pl.DataFrame(rows, schema=schema, orient="row").write_parquet(tmp / "m.parquet")
    t_dense = request.param
    params = IndexParams(k=K, t_base=0.5 if t_dense else 1.0, t_dense=t_dense, fp_bits=64)
    build_index(tmp / "m.parquet", tmp / "idx", params, postings_parquet=True)
    # The "genome": proteins 0-9 whole, the first half of protein 10, reverse strand for 5-9.
    rc = str.maketrans("ACGT", "TGCA")
    parts = [c if i < 5 else c.translate(rc)[::-1] for i, c in enumerate(cds[:10])]
    genome = "NNNN".join([*parts, cds[10][: len(cds[10]) // 2]])
    (tmp / "g.fa").write_text(
        ">chr\n" + "\n".join(genome[i : i + 60] for i in range(0, len(genome), 60)) + "\n"
    )
    build_mask(tmp / "g.fa", tmp / "idx", tmp / "mask")
    return tmp / "idx", tmp / "mask"


def test_decrements_equal_dropping_the_masked_postings(masked: tuple[Path, Path]) -> None:
    idx, mask_dir = masked
    index, mask = Index.load(idx), Mask(mask_dir)
    rows = mask.adjust(index.units.rows(np.arange(index.units.height)))
    postings = pl.read_parquet(idx / "postings.parquet")
    left = postings.filter(~pl.col("hash").is_in(pl.Series(mask.hashes).implode()))
    expected = left.group_by("unit").agg(m=pl.len()).sort("unit")
    got = rows.join(expected, on="unit", how="left").with_columns(pl.col("m").fill_null(0))
    # Units masked whole keep m_g = 1 (they cannot be hit); every other unit loses exactly
    # its masked postings.
    assert (got["m_g"].to_numpy() == np.maximum(got["m"].to_numpy(), 1)).all()
    whole = got.filter(pl.col("m") == 0)["unit"].to_list()
    assert whole == list(range(10))
    host_like = got["host_like"].to_list()
    assert host_like[:10] == [True] * 10
    assert host_like[11:] == [False] * 9
    partly = got.row(10, named=True)
    assert 0 < partly["masked_fraction"] < 1
    hist = np.stack(rows.filter(pl.col("unit") >= 10)["pin_hist"].to_numpy())
    assert (hist.sum(axis=1) == got.filter(pl.col("unit") >= 10)["m_g"].to_numpy()).all()
    if "m_dense" in rows.columns:
        before = index.units.rows(np.array([10]))
        assert rows["m_dense"][10] < before["m_dense"][0]


def test_query_with_the_mask(masked: tuple[Path, Path]) -> None:
    idx, mask_dir = masked
    index = Index.load(idx)
    plain = profile(index, *READS)
    got = profile(index, *READS, mask=Mask(mask_dir))
    assert got.filter(pl.col("unit") < 10).height == 0  # every k-mer of theirs is masked
    assert got.filter(pl.col("unit") < 10).height < plain.filter(pl.col("unit") < 10).height
    # Units the genome does not encode are unchanged, but for the sample-wide gather order
    # and presence prior, which the masked units no longer enter.
    rest = plain.filter(pl.col("unit") > 10).drop("gather_rank", "present_prob")
    assert got.filter(pl.col("unit") > 10).select(rest.columns).equals(rest)


def test_reads_that_hit_nothing(tmp_path: Path) -> None:
    rng = np.random.default_rng(3)
    seqs = ["".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 200)) for _ in range(3)]
    rows = [(i, i, True, q) for i, q in enumerate(seqs)]
    schema = ["protein_id", "cluster_rep", "full_length", "sequence"]
    pl.DataFrame(rows, schema=schema, orient="row").write_parquet(tmp_path / "m.parquet")
    build_index(tmp_path / "m.parquet", tmp_path / "idx", IndexParams(k=K, t_base=1.0))
    got = profile(Index.load(tmp_path / "idx"), *READS, draws=5, with_aai=True)
    assert got.height == 0


def test_cli_mask(masked: tuple[Path, Path], tmp_path: Path) -> None:
    idx, mask_dir = masked
    runner = CliRunner()
    genome = mask_dir.parent / "g.fa"
    result = runner.invoke(app, ["mask", str(genome), str(idx), str(tmp_path / "m")])
    assert result.exit_code == 0, result.output
    assert (tmp_path / "m" / "mask.npy").read_bytes() == (mask_dir / "mask.npy").read_bytes()
    out = tmp_path / "p.tsv"
    args = ["query", str(idx), *map(str, READS), "--out", str(out), "--mask", str(tmp_path / "m")]
    assert runner.invoke(app, args).exit_code == 0
    assert {"masked_fraction", "host_like"} <= set(pl.read_csv(out, separator="\t").columns)
