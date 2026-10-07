"""kfp-prior: carriage shrinkage and the genome-informed presence update."""

import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from test_genomes import DATA, index_dir, records, write_fasta  # noqa: F401  (index_dir: fixture)
from typer.testing import CliRunner

from kfp_prior.cli import app
from kfp_prior.model import Carriage, build_carriage, update, zero_hit_presence
from kmer_functional_profiler.genomes import annotate_genomes


def test_zero_hit_closed_form() -> None:
    # The plan's worked example: q = 0.95, P_G = 1.
    got = zero_hit_presence(np.array([0.95, 0.95, 1 - 1e-6]), np.array([0.5, 5.0, 50.0]))
    assert np.allclose(got[:2], [0.92, 0.11], atol=0.005)
    assert got[2] < 1e-12  # (almost) always carried, 50 hits expected, none seen: absent


def fake_index(tmp: Path, carried: dict[str, set[int]], taxonomy: dict[str, str]) -> Path:
    """A genome index with only what ``build`` reads: genomes.tsv, genome_best, meta.json."""
    tmp.mkdir(parents=True, exist_ok=True)
    names = list(carried)
    pl.DataFrame(
        {"genome": range(len(names)), "name": names, "taxonomy": [taxonomy[n] for n in names]}
    ).write_csv(  # fmt: skip
        tmp / "genomes.tsv", separator="\t"
    )
    rows = [(g, i, u) for g, n in enumerate(names) for i, u in enumerate(sorted(carried[n]))]
    pl.DataFrame(rows, schema=["genome", "protein", "unit"], orient="row").write_parquet(
        tmp / "genome_best.parquet"
    )
    (tmp / "meta.json").write_text(json.dumps({"units": 1000}))
    return tmp


def q_of(
    tmp: Path, name: str, carried: dict[str, set[int]], taxonomy: dict[str, str], alpha: float
) -> dict[int, float]:
    build_carriage(fake_index(tmp / "gi", carried, taxonomy), tmp / "c")
    meta = json.loads((tmp / "c" / "meta.json").read_text())
    meta["alpha"] = dict.fromkeys(meta["alpha"], alpha)  # mechanics; the fit is tested below
    (tmp / "c" / "meta.json").write_text(json.dumps(meta))
    q = Carriage(tmp / "c").q([list(carried).index(name)])
    return dict(zip(q["unit"].to_list(), q["q"].to_list(), strict=True))


def test_shrinkage_towards_one_with_more_genomes_and_towards_genus(tmp_path: Path) -> None:
    small = {f"s{i}": {1, 2} for i in range(2)} | {"o": {3}}
    big = {f"s{i}": {1, 2} for i in range(20)} | {"o": {3}}
    taxonomy = {n: "d__B;f__F;g__G;s__S" if n.startswith("s") else "d__B;f__F;g__G;s__O"
                for n in big}  # fmt: skip
    q2 = q_of(tmp_path / "2", "s0", small, taxonomy, 1.0)
    q20 = q_of(tmp_path / "20", "s0", big, taxonomy, 1.0)
    assert q2[1] < q20[1] and q20[1] > 0.95
    # A one-genome species is pulled half-way (alpha 1) towards its genus, where 1 in 21
    # genomes carries unit 3.
    assert q_of(tmp_path / "o", "o", big, taxonomy, 1.0)[3] < 0.6


def test_held_out_loss_beats_no_shrinkage_and_parent_alone(tmp_path: Path) -> None:
    rng = np.random.default_rng(0)
    carried, taxonomy = {}, {}
    for g in range(6):  # genera with a core, species with their own accessory sets
        core = set(rng.choice(1000, 40, replace=False).tolist())
        for s in range(4):
            accessory = set(rng.choice(1000, 30, replace=False).tolist())
            for i in range(int(rng.integers(2, 6))):
                kept = {u for u in core | accessory if rng.random() < 0.85}
                carried[f"g{g}s{s}_{i}"] = kept
                taxonomy[f"g{g}s{s}_{i}"] = f"d__B;f__F{g % 2};g__G{g};s__S{g}_{s}"
    meta = build_carriage(fake_index(tmp_path / "gi", carried, taxonomy), tmp_path / "c")
    losses = meta["held_out_log_loss"]["species"]  # type: ignore[index]
    assert losses["best"] < losses["no_shrinkage"] and losses["best"] < losses["parent_only"]


@pytest.fixture
def prior_index(tmp_path: Path, index_dir: Path) -> tuple[Path, Path, Path]:  # noqa: F811
    """Species A: a1, a2 (units 0-9), a3 (0-8); species B: b1 (10-19)."""
    seqs = records(DATA / "proteins.faa")
    genomes = {"a1": range(10), "a2": range(10), "a3": range(9), "b1": range(10, 20)}
    lines = ["genome\tpath\ttaxonomy"]
    for g, members in genomes.items():
        write_fasta(tmp_path / f"{g}.faa", [seqs[i] for i in members])
        lines.append(f"{g}\t{g}.faa\td__B;f__F;g__G;s__{g[0]}")
    (tmp_path / "g.tsv").write_text("\n".join(lines) + "\n")
    annotate_genomes(index_dir, tmp_path / "g.tsv", tmp_path / "gi")
    build_carriage(tmp_path / "gi", tmp_path / "c")
    return tmp_path, tmp_path / "gi", tmp_path / "c"


def observed(units: list[int], prob: float = 0.6) -> pl.DataFrame:
    return pl.DataFrame(
        {"unit": units, "hits": [5] * len(units), "present_prob": [prob] * len(units)}
    )


def test_update(prior_index: tuple[Path, Path, Path]) -> None:
    tmp, gi, c = prior_index
    carriage = Carriage(c)
    profile = observed([0, 1, 2, 3, 4, 5, 6, 15])
    # No genomes detected (an empty genomes.tsv reads as strings): every hit unit keeps
    # its present_prob, nothing is imputed.
    none = pl.DataFrame(schema={"genome": pl.Int64, "depth": pl.String, "present": pl.String})
    presence, _ = update(profile, none, gi, carriage)
    assert np.allclose(presence["present_prob_updated"], profile["present_prob"])
    assert presence.height == profile.height
    # Nothing confidently present in the profile: zeros stay zeros, no division by zero.
    genomes = pl.DataFrame({"genome": [0], "depth": [1.0], "present": [1.0]})
    presence, _ = update(observed([0, 1], prob=0.0), genomes, gi, carriage)
    assert presence.filter(pl.col("hits") > 0)["present_prob_updated"].to_list() == [0.0, 0.0]
    # Genome a1 detected. Units 7-9 have no hits: imputed at low depth (unit 9 less, since
    # a3 lacks it), confident absences at high depth.
    for depth, expect in ((1e-4, True), (1.0, False)):
        genomes = pl.DataFrame({"genome": [0], "depth": [depth], "present": [1.0]})
        presence, by_pfam = update(profile, genomes, gi, carriage)
        p = dict(presence.select("unit", "present_prob_updated").iter_rows())
        assert all((p[u] >= 0.5) == expect for u in (7, 8)), p
        assert p[9] < p[7] or not expect
        assert presence.filter(pl.col("unit").is_in([7, 8]))["imputed"].to_list() == [expect] * 2
    # Hit units a1 carries gain; unit 15 (species B, not detected) keeps its value.
    assert p[0] > 0.6 and p[15] == pytest.approx(0.6)
    assert by_pfam is not None and set(by_pfam.columns) >= {"present_prob_observed", "imputed"}
    # Unit 9 shares PF00004 with unit 8 (two units per Pfam in the fixture): no sibling
    # observed; unit 7 shares PF00003 with unit 6, which is observed.
    sib = dict(presence.select("unit", "sibling_hit").iter_rows())
    assert sib[7] and not sib[9]


def test_cli_and_checksum(prior_index: tuple[Path, Path, Path], tmp_path: Path) -> None:
    tmp, gi, c = prior_index
    observed([0, 1]).write_csv(tmp / "p.tsv", separator="\t")
    pl.DataFrame({"genome": [0], "depth": [1e-4], "present": [1.0]}).write_csv(
        tmp / "genomes.tsv", separator="\t"
    )
    args = ["update", str(tmp / "p.tsv"), str(tmp / "genomes.tsv"), str(gi), str(c),
            "--out", str(tmp / "presence.tsv"), "--pfam-out", str(tmp / "pfam.tsv")]  # fmt: skip
    got = CliRunner().invoke(app, args)
    assert got.exit_code == 0, got.output
    assert pl.read_csv(tmp / "presence.tsv", separator="\t")["imputed"].sum() > 0
    other = fake_index(tmp_path / "other", {"x": {1}}, {"x": "s__X"})
    with pytest.raises(ValueError, match="another genome index"):
        update(observed([0]), pl.read_csv(tmp / "genomes.tsv", separator="\t"), other, Carriage(c))
