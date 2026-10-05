"""Genome mode (phase 11): annotation, the genome fit and the function x taxon table."""

from pathlib import Path

import numpy as np
import polars as pl
import pytest
from typer.testing import CliRunner

from kmer_functional_profiler.cli import app
from kmer_functional_profiler.genomes import GenomeIndex, annotate_genomes
from kmer_functional_profiler.index import Index, IndexParams, build_index
from kmer_functional_profiler.query import profile

DATA = Path(__file__).resolve().parents[1] / "data"
K = 7


def records(path: Path) -> list[str]:
    return [r.split("\n", 1)[1].replace("\n", "") for r in path.read_text().split(">")[1:]]


def write_fasta(path: Path, seqs: list[str]) -> Path:
    path.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs)))
    return path


@pytest.fixture(scope="module")
def index_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    tmp = tmp_path_factory.mktemp("genomes")
    rows = [(i, i, True, q) for i, q in enumerate(records(DATA / "proteins.faa"))]
    schema = ["protein_id", "cluster_rep", "full_length", "sequence"]
    pl.DataFrame(rows, schema=schema, orient="row").write_parquet(tmp / "m.parquet")
    pfam = [(i, f"PF{i // 2:05d}") for i in range(20)]  # two units per Pfam
    pl.DataFrame(pfam, schema=["protein_id", "pfam_accession"], orient="row").write_parquet(
        tmp / "pfam.parquet"
    )
    build_index(tmp / "m.parquet", tmp / "idx", IndexParams(k=K, t_base=0.5), tmp / "pfam.parquet")
    return tmp / "idx"


def annotate(tmp: Path, index_dir: Path, genomes: dict[str, list[int]], name: str) -> GenomeIndex:
    """Genome index of genomes made of the fixture proteins listed."""
    seqs = records(DATA / "proteins.faa")
    tsv = ["genome\tpath\ttaxonomy"]
    for g, members in genomes.items():
        write_fasta(tmp / f"{g}.faa", [seqs[i] for i in members])
        tsv.append(f"{g}\t{g}.faa\td__B;g__{g[0]};s__{g}")
    (tmp / f"{name}.tsv").write_text("\n".join(tsv) + "\n")
    annotate_genomes(index_dir, tmp / f"{name}.tsv", tmp / name)
    return GenomeIndex(tmp / name)


def test_genome_of_member_proteins_carries_their_units(tmp_path: Path, index_dir: Path) -> None:
    gi = annotate(tmp_path, index_dir, {"a": [0, 1, 2], "b": [5, 6]}, "gi")
    content = gi.content(np.arange(20))
    carried = content.group_by("genome").agg(pl.col("unit").sort()).sort("genome")
    assert carried["unit"].to_list() == [[0, 1, 2], [5, 6]]
    assert (content["hits"] >= content["kmers"]).all() and (content["kmers"] > 0).all()
    assert gi.genomes["units"].to_list() == [3, 2]
    assert gi.unit_pfam is not None
    gi.check_index(index_dir)
    with pytest.raises(ValueError, match="another index"):
        gi.check_index(tmp_path / "gi")  # any other meta.json
    hits, at = gi.genome_values(np.array([1]))
    assert sorted(hits.tolist()) == sorted(content.filter(pl.col("genome") == 1)["hits"].to_list())
    assert (at == 0).all()


def test_genomes_annotated_together_equal_each_alone(tmp_path: Path, index_dir: Path) -> None:
    together = annotate(tmp_path, index_dir, {"a": [0, 1, 2, 3], "b": [3, 4, 9]}, "both")
    alone = [annotate(tmp_path, index_dir, {g: m}, g) for g, m in (("a", [0, 1, 2, 3]),)]
    alone.append(annotate(tmp_path, index_dir, {"b": [3, 4, 9]}, "b_only"))
    units = np.arange(20)
    for g, single in enumerate(alone):
        expected = single.content(units).drop("genome")
        got = together.content(units).filter(pl.col("genome") == g).drop("genome")
        assert got.equals(expected)


def test_reads_at_depth_hit_units_in_proportion_to_content(tmp_path: Path, index_dir: Path) -> None:
    members = [0, 1, 2, 3, 4]
    gi = annotate(tmp_path, index_dir, {"a": members}, "gi")
    # Error-free 150 nt reads tiled every 3 nt along each CDS (stop codon cut, which would
    # drop the reads over it), padded with N: every k-mer is in (150 - 21) / 3 + 1 = 44
    # reads, so unit hits = 44 x content.
    cds = ["N" * 150 + c[:-3] + "N" * 150 for c in records(DATA / "cds.fna")]
    reads = [c[s : s + 150] for i in members for c in [cds[i]] for s in range(0, len(c) - 149, 3)]
    write_fasta(tmp_path / "r.fa", reads)
    prof = profile(Index.load(index_dir), tmp_path / "r.fa")
    c = gi.content(np.arange(20)).select("unit", c="hits")
    ratio = prof.join(c, on="unit").select(pl.col("hits") / pl.col("c"))["hits"].to_numpy()
    assert len(ratio) == len(members)
    assert np.allclose(ratio, 44, rtol=1e-3)


def test_annotate_genomes_cli(tmp_path: Path, index_dir: Path) -> None:
    write_fasta(tmp_path / "a.faa", records(DATA / "proteins.faa")[:2])
    (tmp_path / "g.tsv").write_text("genome\tpath\na\ta.faa\n")
    got = CliRunner().invoke(app, ["annotate-genomes", str(index_dir), str(tmp_path / "g.tsv"),
                                   str(tmp_path / "out")])  # fmt: skip
    assert got.exit_code == 0, got.output
    assert GenomeIndex(tmp_path / "out").genomes["units"].to_list() == [2]
