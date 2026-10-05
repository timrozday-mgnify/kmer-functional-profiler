"""Genome mode (phase 11): annotation, the genome fit and the function x taxon table."""

import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from typer.testing import CliRunner

from kmer_functional_profiler.cli import app
from kmer_functional_profiler.genomes import GenomeIndex, annotate_genomes, genome_profile
from kmer_functional_profiler.index import Index, IndexParams, build_index
from kmer_functional_profiler.query import em, profile

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
    table, summary, _ = genome_profile(prof, gi, min_units=3)
    assert np.allclose(table["depth"], 44, rtol=1e-3) and np.isclose(table["present"][0], 1.0)
    assert summary["explained_fraction"] > 0.999


def test_annotate_genomes_cli(tmp_path: Path, index_dir: Path) -> None:
    write_fasta(tmp_path / "a.faa", records(DATA / "proteins.faa")[:2])
    (tmp_path / "g.tsv").write_text("genome\tpath\na\ta.faa\n")
    got = CliRunner().invoke(app, ["annotate-genomes", str(index_dir), str(tmp_path / "g.tsv"),
                                   str(tmp_path / "out")])  # fmt: skip
    assert got.exit_code == 0, got.output
    assert GenomeIndex(tmp_path / "out").genomes["units"].to_list() == [2]


def test_weighted_em_with_unit_weights_equals_em() -> None:
    rng = np.random.default_rng(1)
    pairs = pl.DataFrame({"unit": rng.integers(0, 6, 60), "hash": rng.integers(0, 25, 60)}).unique()
    kmers = pairs.with_columns(hits=pl.Series(rng.integers(1, 9, pairs.height)))
    kmers = kmers.with_columns(pl.col("hits").first().over("hash"))  # hits are per k-mer
    m_g = np.full(6, 30)
    for zi in (False, True):
        plain = em(kmers, m_g, zero_inflated=zi)
        weighted = em(kmers.with_columns(weight=pl.lit(1.0)), m_g, zero_inflated=zi)
        assert np.allclose(plain["coverage"], weighted["coverage"])
        assert np.allclose(plain["present"], weighted["present"])


def synthetic(
    gi: GenomeIndex, depth: dict[int, float], units: set[int] | None = None
) -> pl.DataFrame:
    """A profile whose raw unit hits are exactly Σ_G depth_G c_{G,u} (on ``units`` only)."""
    c = gi.content(np.arange(20)).filter(pl.col("genome").is_in(list(depth)))
    if units is not None:
        c = c.filter(pl.col("unit").is_in(list(units)))
    weights = pl.DataFrame({"genome": list(depth), "d": list(depth.values())})
    return (
        c.with_columns(pl.col("genome").cast(pl.Int64))
        .join(weights, on="genome")
        .group_by("unit")
        .agg(hits=(pl.col("d") * pl.col("hits")).sum().round())
        .sort("unit")
    )


def test_disjoint_genomes_recover_their_depths(tmp_path: Path, index_dir: Path) -> None:
    genomes = {"a": list(range(10)), "b": list(range(10, 20)), "x": [0, 1, 2, 10, 11]}
    gi = annotate(tmp_path, index_dir, genomes, "gi")
    table, summary, _ = genome_profile(synthetic(gi, {0: 3.0, 1: 7.0}), gi, min_units=5)
    assert table["name"].to_list() == ["b", "a"]  # x: explained away by gather
    assert np.allclose(table["depth"], [7.0, 3.0], rtol=0.01)
    assert np.allclose(table["present"], 1.0, atol=0.01)
    assert summary["explained_fraction"] > 0.99
    assert np.isclose(summary["genome_equivalents"], 10.0, rtol=0.01)


def test_relative_missing_content_gets_present_below_one(tmp_path: Path, index_dir: Path) -> None:
    """The sample's strain is a relative of genome a lacking 3 of its 10 units."""
    gi = annotate(tmp_path, index_dir, {"a": list(range(10))}, "gi")
    table, _, _ = genome_profile(synthetic(gi, {0: 5.0}, set(range(7))), gi, min_units=5)
    assert table["name"].to_list() == ["a"]
    assert 0.6 < table["present"][0] < 0.8
    assert abs(table["depth"][0] - 5.0) / 5.0 < 0.2


def test_unexplained_hits_and_cli(tmp_path: Path, index_dir: Path) -> None:
    gi = annotate(tmp_path, index_dir, {"a": list(range(10))}, "gi")
    prof = pl.concat(
        [
            synthetic(gi, {0: 2.0}),
            pl.DataFrame({"unit": [15], "hits": [100.0]}).cast({"unit": pl.UInt32}),
        ]
    )
    prof.write_csv(tmp_path / "p.tsv", separator="\t")
    args = ["genomes", str(tmp_path / "p.tsv"), str(tmp_path / "gi"), str(tmp_path / "g.tsv")]
    args += ["--summary", str(tmp_path / "s.json"), "--index", str(index_dir), "--min-units", "5"]
    got = CliRunner().invoke(app, args)
    assert got.exit_code == 0, got.output
    summary = json.loads((tmp_path / "s.json").read_text())
    hits = summary["hits"]
    assert np.isclose(summary["explained_fraction"], (hits - 100) / hits, rtol=0.01)
    bad = CliRunner().invoke(app, [*args[:6], "--index", str(tmp_path / "gi")])
    assert bad.exit_code != 0


def stratified_profile(
    gi: GenomeIndex, depth: dict[int, float], extra: dict[int, float]
) -> pl.DataFrame:
    """:func:`synthetic` with ``hits_em`` = ``hits`` and ``extra`` hits on other units."""
    prof = pl.concat(
        [
            synthetic(gi, depth),
            pl.DataFrame({"unit": list(extra), "hits": list(extra.values())}).cast(
                {"unit": pl.UInt32}
            ),
        ]
    )
    return prof.with_columns(hits_em=pl.col("hits"), name=pl.col("unit").cast(pl.String))


def test_function_taxon_table(tmp_path: Path, index_dir: Path) -> None:
    # a: units 0-8, b: units 10-18 and unit 0; units 9 and 19 carried by neither.
    genomes = {"a": list(range(9)), "b": [0, *range(10, 19)]}
    gi = annotate(tmp_path, index_dir, genomes, "gi")
    c = gi.content(np.arange(20))
    prof = stratified_profile(gi, {0: 1.0, 1: 3.0}, {19: 50.0})
    table, _, ft = genome_profile(prof, gi, min_units=5)
    assert ft is not None
    assert np.allclose(table.sort("genome")["depth"], [1.0, 3.0], rtol=0.01)
    assert table["group"].n_unique() == 2
    assert ((table["depth_lo"] <= table["depth"]) & (table["depth"] <= table["depth_hi"])).all()

    def pf(u: int) -> str:  # the fixture's labels: two units per Pfam
        return f"PF{u // 2:05d}"

    sp = ft.filter(pl.col("rank") == "species")
    # A unit of a alone (unit 2, PF00001 with unit 3, also a's) goes to a whole.
    one = sp.filter(pl.col("function") == pf(2))
    assert one["taxon"].to_list() == ["s__a"]
    # Unit 0 (shared, PF00000 with a's unit 1) splits 1:3 between a and b.
    c0 = c.filter(pl.col("unit") == 0)["hits"][0]
    c1 = c.filter(pl.col("unit") == 1)["hits"][0]
    got = dict(sp.filter(pl.col("function") == pf(0)).select("taxon", "hits_em").iter_rows())
    assert np.isclose(got["s__a"], 1.0 * c0 + c1, rtol=0.01)
    assert np.isclose(got["s__b"], 3.0 * c0, rtol=0.01)
    # Unit 19 (PF00009 with unit 18, b's) is carried by no genome: unclassified.
    got = dict(sp.filter(pl.col("function") == pf(19)).select("taxon", "hits_em").iter_rows())
    assert np.isclose(got["unclassified"], 50.0)
    # Every rank's rows sum to the total.
    totals = ft.filter(pl.col("rank") == "total").select("function", total="hits_em")
    sums = (
        ft.filter(pl.col("rank") != "total")
        .group_by("function", "rank")
        .agg(pl.col("hits_em").sum())
    )
    check = sums.join(totals, on="function")
    assert set(sums["rank"]) == {"family", "genus", "species", "genome"}
    assert np.allclose(check["hits_em"], check["total"])


def test_identical_genomes_form_a_group_reported_at_their_common_rank(
    tmp_path: Path, index_dir: Path
) -> None:
    gi = annotate(tmp_path, index_dir, {"a1": list(range(10)), "a2": list(range(10))}, "gi")
    prof = stratified_profile(gi, {0: 4.0}, {})
    table, summary, ft = genome_profile(prof, gi, min_units=5)
    assert ft is not None and summary["ambiguity_groups"] == 1
    assert np.allclose(table["depth"], 4.0, rtol=0.01)
    low = ft.filter(pl.col("rank").is_in(["species", "genome"]))["taxon"].unique()
    assert low.to_list() == ["g__a"]  # genus shared; species and genome differ
    assert table["ambiguous"].to_list() == [1]  # a2, explained away by gather
