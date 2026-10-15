"""The read cross-reference join (kfp_genomes.xref) on hand-made tables."""

import polars as pl
import pytest

from kfp_genomes.xref import HITS_COLUMNS, cross_reference


def _assignments(rows: list[tuple[str, int, int, str, float, bool]]) -> pl.DataFrame:
    return pl.DataFrame(
        rows,
        schema={
            "read_name": pl.String,
            "mate": pl.UInt8,
            "unit": pl.UInt32,
            "role": pl.String,
            "posterior": pl.Float64,
            "kept": pl.Boolean,
        },
        orient="row",
    )


def test_join_splits_pairs_and_keeps_unclassified() -> None:
    assignments = _assignments(
        [
            # p1: both mates on unit 0 (a decoy runner-up is renormalised away)
            ("p1/1", 0, 0, "target", 0.9, True),
            ("p1/1", 0, 0, "decoy", 0.1, True),
            ("p1/2", 1, 0, "target", 1.0, True),
            # p2: mates disagree, units 0 and 1
            ("p2/1", 0, 0, "target", 1.0, True),
            ("p2/2", 1, 1, "target", 1.0, True),
            # p3: unit 1, no genome hit -> unclassified
            ("p3/1", 0, 1, "target", 1.0, True),
            # p4: best hit a neighbour: dropped
            ("p4/1", 0, 1, "neighbour", 1.0, False),
        ]
    )
    refs = pl.DataFrame(
        {
            "unit": [0, 0, 1],
            "role": ["target"] * 3,
            "length": [100, 200, 50],
            "rep": [True, False, True],
        },
        schema_overrides={"unit": pl.UInt32},
    )
    hits = pl.DataFrame(
        [
            # p1: both mates on gA; one mate on gB only -> gA wins on mates
            ("p1", 65, "gA|0", 0, -2, 1, 100),
            ("p1", 129, "gA|0", 0, -3, 0, 100),
            ("p1", 65, "gB|0", 0, 0, 0, 100),
            # p2: tied between gA and gB -> split; a low-identity hit on gC is dropped
            ("p2", 65, "gA|0", 0, 0, 0, 100),
            ("p2", 65, "gB|0", 0, 0, 0, 100),
            ("p2", 65, "gC|0", 0, 10, 20, 100),
            # p5: genome only (no function)
            ("p5", 65, "gB|0", 0, 0, 0, 100),
            ("p6", 4, "*", 0, 0, 0, 0),  # unaligned
        ],
        schema=HITS_COLUMNS,
        orient="row",
    )
    contigs = pl.DataFrame(
        {
            "contig": ["gA|0", "gB|0", "gC|0"],
            "genome": ["gA", "gB", "gC"],
            "length": [1000, 2000, 500],
        }
    )
    taxonomy = pl.DataFrame(
        {"genome": ["gA", "gB", "gC"], "taxonomy": ["g__X;s__X a", "g__X;s__X b", "g__Y;s__Y c"]}
    )
    fg, genomes = cross_reference(assignments, refs, False, hits, contigs, taxonomy)

    g = fg.filter(pl.col("rank") == "genome")
    reads = {(r["unit"], r["taxon"]): r["reads"] for r in g.iter_rows(named=True)}
    assert reads == pytest.approx(
        {
            (0, "gA"): 1.25,
            (0, "gB"): 0.25,
            (1, "gA"): 0.25,
            (1, "gB"): 0.25,
            (1, "unclassified"): 1.0,
        }
    )
    # rows of a unit and rank sum to the unit's pairs; share sums to 1
    totals = fg.group_by("unit", "rank").agg(pl.col("reads", "share").sum()).sort("unit", "rank")
    assert totals["reads"].to_list() == pytest.approx([1.5] * 3 + [1.5] * 3)
    assert totals["share"].to_list() == pytest.approx([1.0] * 6)
    # species gA, gB both in genus X
    genus = fg.filter(pl.col("rank") == "genus", pl.col("unit") == 0)
    assert genus["taxon"].to_list() == ["X"] and genus["reads"][0] == pytest.approx(1.5)

    # genome depth: p1 200 bases on gA, p2 half of 100 each, p5 100 on gB
    by = {r["taxon"]: r for r in genomes.filter(pl.col("rank") == "genome").iter_rows(named=True)}
    assert by["gA"]["reads"] == pytest.approx(1.5) and by["gB"]["reads"] == pytest.approx(1.5)
    assert by["gA"]["depth"] == pytest.approx(250 / 1000)
    assert by["gB"]["depth"] == pytest.approx(150 / 2000)
    assert by["gB"]["function_reads"] == pytest.approx(0.5) and by["gC"]["reads"] == 0
    # coverage uses the representative's length (100 aa = 300 bases); copies = coverage / depth
    row = g.filter(pl.col("unit") == 0, pl.col("taxon") == "gA").row(0, named=True)
    bases_per_pair = (200 + 50 + 50 + 100) / 3
    assert row["coverage"] == pytest.approx(1.25 * bases_per_pair / 300)
    assert row["copies"] == pytest.approx(row["coverage"] / 0.25)
    assert g.filter(pl.col("taxon") == "unclassified")["copies"].is_null().all()
