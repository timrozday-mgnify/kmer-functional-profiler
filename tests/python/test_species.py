"""Species model (phase 11, steps 7-9): the species index, the fit and its outputs."""

# ruff: noqa: F811  (index_dir: the fixture imported from test_genomes)

import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from test_genomes import DATA, index_dir, records, write_fasta  # noqa: F401  (index_dir: fixture)
from typer.testing import CliRunner

from kmer_functional_profiler.cli import app
from kmer_functional_profiler.genomes import GenomeIndex, annotate_genomes
from kmer_functional_profiler.species import (
    RANKS,
    SpeciesIndex,
    _counts,
    fit_alpha,
    lineages,
    prevalence,
    species_index_from_catalogue,
    species_index_from_genomes,
    translate,
)

# --- Shrinkage (moved from kfp-prior) ---


def q_of(carried: dict[str, set[int]], taxonomy: dict[str, str], alpha: float) -> pl.DataFrame:
    """Per species and unit, q at a fixed α (the α fit is tested below)."""
    names = list(carried)
    genomes = pl.DataFrame(
        {"genome": range(len(names)), "name": names, "taxonomy": [taxonomy[n] for n in names]}
    ).cast({"genome": pl.UInt32})
    rows = [(g, u, 1.0) for g, n in enumerate(names) for u in sorted(carried[n])]
    rows_df = pl.DataFrame(rows, schema=["genome", "unit", "c"], orient="row")
    c = rows_df.cast({"genome": pl.UInt32, "unit": pl.UInt32})
    lineage = lineages(genomes)
    n, sizes = _counts(c.select("genome", "unit"), lineage)
    return prevalence(c, lineage, n, sizes, dict.fromkeys(RANKS, alpha))


def q_at(table: pl.DataFrame, species: str, unit: int) -> float:
    row = table.filter((pl.col("species") == species) & (pl.col("unit") == unit))
    return float(row["q"][0]) if row.height else 0.0


def test_shrinkage_towards_one_with_more_genomes_and_towards_genus() -> None:
    small = {f"s{i}": {1, 2} for i in range(2)} | {"o": {3}}
    big = {f"s{i}": {1, 2} for i in range(20)} | {"o": {3}}
    taxonomy = {n: "d__B;f__F;g__G;s__S" if n.startswith("s") else "d__B;f__F;g__G;s__O"
                for n in big}  # fmt: skip
    q2, q20 = q_of(small, taxonomy, 1.0), q_of(big, taxonomy, 1.0)
    assert q_at(q2, "s__S", 1) < q_at(q20, "s__S", 1) and q_at(q20, "s__S", 1) > 0.95
    # A one-genome species is pulled half-way (alpha 1) towards its genus, where 1 in 21
    # genomes carries unit 3; and it gets a (small) prior on its genus' units 1 and 2.
    assert q_at(q20, "s__O", 3) < 0.6
    assert 0.4 < q_at(q20, "s__O", 1) < 0.5


def test_held_out_loss_beats_no_shrinkage_and_parent_alone() -> None:
    rng = np.random.default_rng(0)
    rows, taxonomy = [], {}
    g_id = 0
    for g in range(6):  # genera with a core, species with their own accessory sets
        core = set(rng.choice(1000, 40, replace=False).tolist())
        for s in range(4):
            accessory = set(rng.choice(1000, 30, replace=False).tolist())
            for _ in range(int(rng.integers(2, 6))):
                rows += [(g_id, u) for u in core | accessory if rng.random() < 0.85]
                taxonomy[g_id] = f"d__B;f__F{g % 2};g__G{g};s__S{g}_{s}"
                g_id += 1
    genomes = pl.DataFrame(
        {"genome": list(taxonomy), "name": [str(g) for g in taxonomy],
         "taxonomy": list(taxonomy.values())}
    ).cast({"genome": pl.UInt32})  # fmt: skip
    carried = pl.DataFrame(rows, schema=["genome", "unit"], orient="row").cast(
        {"genome": pl.UInt32}
    )
    _, losses = fit_alpha(carried, lineages(genomes), "species")
    assert losses["best"] < losses["no_shrinkage"] and losses["best"] < losses["parent_only"]


# --- Species index from a catalogue and from a genome set ---

# Species A (rep MGYG000000001): three genomes over families 0-9; genome 3 (50% complete)
# lacks families 8 and 9, genome 2 lacks family 9. Species B (MGYG000000004): one genome,
# proteins 10-19, no pangenome.
A_GENOMES = {"MGYG000000001": range(10), "MGYG000000002": range(9), "MGYG000000003": range(8)}
COMPLETENESS = {"MGYG000000001": 100, "MGYG000000002": 100, "MGYG000000003": 50,
                "MGYG000000004": 100}  # fmt: skip
LINEAGE = "d__B;p__P;c__C;o__O;f__F;g__G;s__"


def cds_records() -> list[str]:
    return records(DATA / "cds.fna")


def fake_catalogue(tmp: Path, completeness: dict[str, int] = COMPLETENESS) -> Path:
    root = tmp / "catalogue"
    a = root / "species_catalogue" / "MGYG0000000" / "MGYG000000001" / "pan-genome"
    a.mkdir(parents=True)
    cds = cds_records()
    (a / "pan-genome.fna").write_text("".join(f">fam{i}\n{cds[i]}\n" for i in range(10)))
    header = "Gene\t" + "\t".join(A_GENOMES)
    rows = [f"fam{i}\t" + "\t".join(str(int(i in m)) for m in A_GENOMES.values())
            for i in range(10)]  # fmt: skip
    (a / "gene_presence_absence.Rtab").write_text("\n".join([header, *rows]) + "\n")
    b = root / "species_catalogue" / "MGYG0000000" / "MGYG000000004" / "genome"
    b.mkdir(parents=True)
    write_fasta(b / "MGYG000000004.faa", records(DATA / "proteins.faa")[10:20])
    meta = ["Genome\tCompleteness\tSpecies_rep\tLineage"]
    for g in A_GENOMES:
        meta.append(f"{g}\t{completeness[g]}\tMGYG000000001\t{LINEAGE}A")
    meta.append(f"MGYG000000004\t{completeness['MGYG000000004']}\tMGYG000000004\t{LINEAGE}B")
    (root / "genomes-all_metadata.tsv").write_text("\n".join(meta) + "\n")
    return root


def pairs(si: SpeciesIndex) -> pl.DataFrame:
    names = si.species.select(pl.col("species").cast(pl.UInt32), "name")
    return (
        si.by_species(np.arange(si.species.height))
        .join(names, on="species")
        .select("name", "unit", "q", "e")
        .sort("name", "unit")
    )


def test_translate_matches_the_proteins() -> None:
    proteins = records(DATA / "proteins.faa")
    assert [translate(c.encode()).decode() for c in cds_records()] == proteins


def test_species_index_from_catalogue(tmp_path: Path, index_dir: Path) -> None:
    species_index_from_catalogue(index_dir, fake_catalogue(tmp_path), tmp_path / "si")
    si = SpeciesIndex(tmp_path / "si")
    si.check_index(index_dir)
    assert si.species["name"].to_list() == ["A", "B"]
    assert si.species["genomes"].to_list() == [3, 1]
    assert si.species["completeness"].to_list() == [2.5, 1.0]
    p = pairs(si)
    a = p.filter(pl.col("name") == "A")
    # Units 0-7: in every genome, 3 / 2.5 > 1, capped. Unit 9: in one genome of 2.5.
    assert np.allclose(a.filter(pl.col("unit") < 8)["q"], 1.0)
    assert (
        q_at(a.rename({"name": "species"}), "A", 9)
        < q_at(a.rename({"name": "species"}), "A", 8)
        < 1
    )
    # e: every carrier's content is its family's, the same protein: the genome-set index's c.
    gi = GenomeIndex(annotate_genomes_for(tmp_path, index_dir, {"x": list(range(20))}))
    c = dict(gi.content(np.arange(20)).select("unit", "hits").iter_rows())
    assert all(np.isclose(e, c[u]) for u, e in p.select("unit", "e").iter_rows())
    # B: one genome, shrunk towards genus G, and given A's units at a small prior.
    b = p.filter(pl.col("name") == "B")
    assert (b.filter(pl.col("unit") >= 10)["q"] < 1).all()
    assert b.filter(pl.col("unit") < 10).height > 0
    assert (b.filter(pl.col("unit") < 10)["q"] < 0.6).all()
    # unit-major view: the same pairs
    by_unit = si.by_unit(np.arange(20)).sort("species", "unit")
    assert by_unit.height == p.height
    assert np.allclose(by_unit["q"], si.by_species(np.arange(2)).sort("species", "unit")["q"])


def test_completeness_raises_prevalence_of_units_incomplete_genomes_miss(
    tmp_path: Path, index_dir: Path
) -> None:
    complete = dict.fromkeys(COMPLETENESS, 100)
    species_index_from_catalogue(index_dir, fake_catalogue(tmp_path), tmp_path / "w")
    species_index_from_catalogue(
        index_dir, fake_catalogue(tmp_path / "c", complete), tmp_path / "u"
    )
    weighted, unweighted = (pairs(SpeciesIndex(tmp_path / x)).rename({"name": "species"})
                            for x in ("w", "u"))  # fmt: skip
    assert q_at(weighted, "A", 8) > q_at(unweighted, "A", 8)


def annotate_genomes_for(tmp: Path, index_dir: Path, genomes: dict[str, list[int]]) -> Path:
    seqs = records(DATA / "proteins.faa")
    lines = ["genome\tpath\ttaxonomy\tcompleteness"]
    for g, members in genomes.items():
        write_fasta(tmp / f"{g}.faa", [seqs[i] for i in members])
        species = "B" if g == "MGYG000000004" else "A"
        lines.append(f"{g}\t{g}.faa\t{LINEAGE}{species}\t{COMPLETENESS.get(g, 100)}")
    (tmp / "g.tsv").write_text("\n".join(lines) + "\n")
    annotate_genomes(index_dir, tmp / "g.tsv", tmp / "gi")
    return tmp / "gi"


def test_genome_set_path_equals_catalogue_path(tmp_path: Path, index_dir: Path) -> None:
    species_index_from_catalogue(index_dir, fake_catalogue(tmp_path), tmp_path / "cat")
    genomes = {g: list(m) for g, m in A_GENOMES.items()} | {"MGYG000000004": list(range(10, 20))}
    gi = annotate_genomes_for(tmp_path, index_dir, genomes)
    species_index_from_genomes(index_dir, gi, tmp_path / "set")
    cat, gset = pairs(SpeciesIndex(tmp_path / "cat")), pairs(SpeciesIndex(tmp_path / "set"))
    assert cat.select("name", "unit").equals(gset.select("name", "unit"))
    assert np.allclose(cat["q"], gset["q"]) and np.allclose(cat["e"], gset["e"])
    with pytest.raises(ValueError, match="another index"):
        species_index_from_genomes(tmp_path / "set", gi, tmp_path / "bad")


def test_exclude_and_cli(tmp_path: Path, index_dir: Path) -> None:
    (tmp_path / "ex.txt").write_text("MGYG000000002\n")
    args = ["species-index", str(index_dir), str(tmp_path / "si"),
            "--catalogue", str(fake_catalogue(tmp_path)),
            "--exclude", str(tmp_path / "ex.txt")]  # fmt: skip
    got = CliRunner().invoke(app, args)
    assert got.exit_code == 0, got.output
    si = SpeciesIndex(tmp_path / "si")
    assert si.species["genomes"].to_list() == [2, 1] and si.meta["excluded"] == 1
    held = pl.read_parquet(tmp_path / "si" / "held_out.parquet")
    assert held["name"].unique().to_list() == ["MGYG000000002"]
    assert held["unit"].to_list() == list(range(9))
    assert json.loads((tmp_path / "si" / "meta.json").read_text())["source"] == "catalogue"
    bad = CliRunner().invoke(app, ["species-index", str(index_dir), str(tmp_path / "x")])
    assert bad.exit_code != 0
