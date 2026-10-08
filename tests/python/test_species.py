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
    fit_species,
    lineages,
    prevalence,
    species_index_from_catalogue,
    species_index_from_genomes,
    species_profile,
    translate,
    write_species_index,
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


# --- The fit and its outputs ---


def synthetic_index(
    tmp: Path, species: dict[str, dict[int, float]], genomes: int = 10, n_units: int = 200
) -> SpeciesIndex:
    """A species index of ``genomes`` genomes per species, each carrying unit u with
    probability species[s][u] (content 2 + u % 7), units 0..max: no profiler index needed."""
    rng = np.random.default_rng(0)
    (tmp / "idx").mkdir(parents=True, exist_ok=True)
    (tmp / "idx" / "meta.json").write_text("{}")
    rows, gen = [], []
    for name, carriage in species.items():
        for _ in range(genomes):
            g = len(gen)
            gen.append((g, f"{name}{g}", f"d__B;f__F;g__G{name};s__{name}"))
            rows += [(g, u, 2.0 + u % 7) for u, p in carriage.items() if rng.random() < p]
    carried = pl.DataFrame(rows, schema=["genome", "unit", "c"], orient="row").cast(
        {"genome": pl.UInt32, "unit": pl.UInt32}
    )
    table = pl.DataFrame(gen, schema=["genome", "name", "taxonomy"], orient="row").cast(
        {"genome": pl.UInt32}
    )
    write_species_index(tmp / "si", tmp / "idx", n_units, carried, table, frozenset(), {})
    return SpeciesIndex(tmp / "si")


def exact(si: SpeciesIndex, depth: dict[str, float], units: dict[str, set[int]]) -> pl.DataFrame:
    """A profile whose hits are exactly Σ_s depth_s e_{s,u} over the units each strain has."""
    ids = dict(zip(si.species["name"], si.species["species"], strict=True))
    hits: dict[int, float] = {}
    for name, d in depth.items():
        pairs = si.by_species(np.array([ids[name]]))
        for u, e in pairs.select("unit", "e").iter_rows():
            if u in units[name]:
                hits[u] = hits.get(u, 0.0) + d * e
    return pl.DataFrame({"unit": list(hits), "hits": list(hits.values())}).cast({"unit": pl.UInt32})


def fitted(si: SpeciesIndex, prof: pl.DataFrame) -> tuple[dict[str, dict], pl.DataFrame]:
    species, pairs, _ = fit_species(prof, si)
    names = si.species.select(pl.col("species").cast(pl.UInt32), "name")
    by_name = {row["name"]: row for row in species.join(names, on="species").iter_rows(named=True)}
    return by_name, pairs.join(names, on="species")


CORE_50 = dict.fromkeys(range(50), 1.0)


def test_one_core_species_recovers_its_depth(tmp_path: Path) -> None:
    si = synthetic_index(tmp_path, {"A": CORE_50})
    by_name, pairs = fitted(si, exact(si, {"A": 2.0}, {"A": set(range(50))}))
    assert by_name["A"]["present_prob"] > 0.999
    assert by_name["A"]["depth"] == pytest.approx(2.0, rel=0.01)
    assert (pairs["carriage_prob"] > 0.999).all()


def test_accessory_units_follow_depth_and_hits(tmp_path: Path) -> None:
    accessory = dict.fromkeys(range(50, 70), 0.5)
    si = synthetic_index(tmp_path, {"A": CORE_50 | accessory})
    strain = set(range(60))  # half the accessory units, as their prevalence has it
    for depth in (5.0, 0.01):
        prof = exact(si, {"A": depth}, {"A": strain})
        by_name, pairs = fitted(si, prof)
        assert by_name["A"]["depth"] == pytest.approx(depth, rel=0.05)
        r = dict(pairs.select("unit", "carriage_prob").iter_rows())
        q = dict(pairs.select("unit", "q").iter_rows())
        if depth > 1:  # zero hits where ~20 were due: absent; hit: carried
            assert all(r[u] < 1e-3 for u in range(60, 70)) and all(r[u] > 0.99 for u in strain)
        else:  # ~0.05 hits due: zero hits say little, carriage stays near prevalence
            assert all(0.9 * q[u] < r[u] <= q[u] for u in range(60, 70)), r
    # At zero hits and no other source, carriage is kfp-prior's closed form (at the
    # posterior mean depth).
    row = pairs.filter(pl.col("unit") == 60)
    q, e, a = row["q"][0], row["e"][0], by_name["A"]
    kept = q * np.exp(-a["shape"] / a["rate"] * e)
    assert row["carriage_prob"][0] == pytest.approx(kept / (1 - q + kept), rel=1e-3)


def test_background_and_a_species_whose_core_is_missing(tmp_path: Path) -> None:
    si = synthetic_index(tmp_path, {"A": CORE_50, "B": dict.fromkeys(range(50, 100), 1.0)})
    prof = exact(si, {"A": 1.0}, {"A": set(range(50))})
    # an organism outside the set hits 15 of B's 50 core units, and units nobody carries
    other = pl.DataFrame({"unit": list(range(50, 65)), "hits": [8.0] * 15})
    by_name, _ = fitted(si, pl.concat([prof, other.cast({"unit": pl.UInt32})]))
    assert by_name["A"]["present_prob"] > 0.999
    assert by_name["A"]["depth"] == pytest.approx(1.0, rel=0.02)
    assert by_name["B"]["present_prob"] < 1e-3


def test_two_species_sharing_core_units(tmp_path: Path) -> None:
    shared = dict.fromkeys(range(30), 1.0)
    si = synthetic_index(
        tmp_path,
        {"A": shared | dict.fromkeys(range(30, 80), 1.0),
         "B": shared | dict.fromkeys(range(80, 130), 1.0)},
    )  # fmt: skip
    units = {"A": set(range(80)), "B": set(range(30)) | set(range(80, 130))}
    by_name, pairs = fitted(si, exact(si, {"A": 1.0, "B": 3.0}, units))
    assert by_name["A"]["depth"] == pytest.approx(1.0, rel=0.02)
    assert by_name["B"]["depth"] == pytest.approx(3.0, rel=0.02)
    # a shared unit's hits split 1:3
    split = dict(pairs.filter(pl.col("unit") == 0).select("name", "hits_assigned").iter_rows())
    assert split["B"] / split["A"] == pytest.approx(3.0, rel=0.02)


def test_profile_and_outputs(tmp_path: Path) -> None:
    accessory = dict.fromkeys(range(50, 70), 0.5)
    si = synthetic_index(tmp_path, {"A": CORE_50 | accessory})
    # A at low depth: misses part of its core; a unit nobody carries has hits.
    rng = np.random.default_rng(1)
    prof = exact(si, {"A": 0.15}, {"A": set(range(50))}).with_columns(
        hits=pl.Series(rng.poisson(exact(si, {"A": 0.15}, {"A": set(range(50))})["hits"]))
    )
    prof = pl.concat([prof, pl.DataFrame({"unit": [70], "hits": [20]})], how="vertical_relaxed")
    prof = (
        prof.cast({"hits": pl.Float64})
        .with_columns(
            unit=pl.col("unit").cast(pl.UInt32), present_prob=pl.lit(0.9), hits_em=pl.col("hits")
        )
        .filter(pl.col("hits") > 0)
    )
    out = species_profile(prof, si, min_units=5)
    species = out["species"]
    assert isinstance(species, pl.DataFrame) and species["name"].to_list() == ["A"]
    assert species["depth_lo"][0] < 0.15 < species["depth_hi"][0]
    presence = out["presence"]
    assert isinstance(presence, pl.DataFrame)
    p = dict(presence.select("unit", "present_prob_updated").iter_rows())
    missed = set(range(50)) - set(prof["unit"].to_list())
    assert missed and all(p[u] > 0.5 for u in missed)  # core units imputed
    assert p[70] == pytest.approx(0.9)  # no species carries it: unchanged
    assert presence.filter(pl.col("unit").is_in(list(missed)))["imputed"].all()
    table = out["function_species"]
    assert isinstance(table, pl.DataFrame)
    totals = table.filter(pl.col("rank") == "total").select("function", total="hits_em")
    sums = (
        table.filter(pl.col("rank") != "total")
        .group_by("function", "rank")
        .agg(pl.col("hits_em").sum())
    )
    check = sums.join(totals, on="function")
    assert set(sums["rank"]) == {"family", "genus", "species"}
    assert np.allclose(check["hits_em"], check["total"])
    assert set(table.filter(pl.col("function") == "70")["taxon"]) == {"", "unclassified"}
    sp = table.filter((pl.col("rank") == "species") & (pl.col("taxon") == "A"))
    assert sp.height == prof.filter(pl.col("unit") < 70).height  # every hit unit of A


def test_no_species_leaves_presence_unchanged(tmp_path: Path) -> None:
    si = synthetic_index(tmp_path, {"A": CORE_50})
    prof = pl.DataFrame({"unit": [3, 7], "hits": [5.0, 2.0], "present_prob": [0.4, 0.8]}).cast(
        {"unit": pl.UInt32}
    )
    out = species_profile(prof, si)
    presence = out["presence"]
    assert isinstance(presence, pl.DataFrame)
    assert out["species"].height == 0  # type: ignore[union-attr]
    hit = presence.filter(pl.col("hits") > 0).sort("unit")
    assert hit["present_prob_updated"].to_list() == pytest.approx([0.4, 0.8])
    assert not presence["imputed"].any()


def test_species_cli(tmp_path: Path, index_dir: Path) -> None:
    genomes = {g: list(m) for g, m in A_GENOMES.items()} | {"MGYG000000004": list(range(10, 20))}
    gi = annotate_genomes_for(tmp_path, index_dir, genomes)
    species_index_from_genomes(index_dir, gi, tmp_path / "si")
    si = SpeciesIndex(tmp_path / "si")
    prof = exact(si, {"A": 3.0}, {"A": set(range(10))}).with_columns(
        hits_em=pl.col("hits"), present_prob=pl.lit(1.0)
    )
    prof.write_csv(tmp_path / "p.tsv", separator="\t")
    out = tmp_path / "o"
    out.mkdir()
    args = ["species", str(tmp_path / "p.tsv"), str(tmp_path / "si"), str(out / "s.tsv"),
            "--units", str(out / "u.tsv"), "--presence", str(out / "pr.tsv"),
            "--pfam-presence", str(out / "pf.tsv"), "--function-taxon", str(out / "ft.tsv"),
            "--summary", str(out / "sum.json"), "--index", str(index_dir),
            "--min-units", "5"]  # fmt: skip
    got = CliRunner().invoke(app, args)
    assert got.exit_code == 0, got.output
    species = pl.read_csv(out / "s.tsv", separator="\t")
    assert species["name"].to_list() == ["A"]
    assert species["depth"][0] == pytest.approx(3.0, rel=0.02)
    ft = pl.read_csv(out / "ft.tsv", separator="\t")
    assert ft["function"].str.starts_with("PF").all()  # the index's Pfam labels
    assert json.loads((out / "sum.json").read_text())["explained_fraction"] > 0.99
    assert pl.read_csv(out / "pf.tsv", separator="\t").height > 0
    bad = CliRunner().invoke(app, [*args[:4], "--index", str(tmp_path / "si")])
    assert bad.exit_code != 0


def test_identical_species_are_explained_away(tmp_path: Path) -> None:
    si = synthetic_index(tmp_path, {"A": CORE_50, "B": CORE_50})
    species, _, report = fit_species(exact(si, {"A": 3.0}, {"A": set(range(50))}), si)
    assert report["screened_species"] == 2 and species.height == 1  # gather keeps one
    assert species["depth"][0] == pytest.approx(3.0, rel=0.01)


def test_overdispersed_hits_keep_the_species_and_its_core(tmp_path: Path) -> None:
    """Gamma-Poisson hits (variance ~6x the mean): Poisson evidence alone would switch
    carriage off on the units far from λe; tempering by φ keeps them."""
    si = synthetic_index(tmp_path, {"A": CORE_50})
    mean = exact(si, {"A": 4.0}, {"A": set(range(50))})
    rng = np.random.default_rng(3)
    hits = rng.poisson(rng.gamma(0.8, mean["hits"].to_numpy() / 0.8))
    prof = mean.with_columns(hits=pl.Series(hits, dtype=pl.Float64)).filter(pl.col("hits") > 0)
    species, pairs, report = fit_species(prof, si)
    assert report["fit_phi"] > 3
    assert species["present_prob"][0] > 0.99
    assert species["depth"][0] == pytest.approx(4.0, rel=0.3)
    hit = pairs.filter(pl.col("hits") > 0)
    assert (hit["carriage_prob"] > 0.9).all()
