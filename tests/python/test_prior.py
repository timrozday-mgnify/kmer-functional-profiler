"""kfp-prior: the genome-informed presence update (shrinkage: test_species.py)."""

import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from test_genomes import DATA, index_dir, records, write_fasta  # noqa: F401  (index_dir: fixture)
from typer.testing import CliRunner

from kfp_prior.cli import app
from kfp_prior.model import Carriage, update, zero_hit_presence
from kmer_functional_profiler.genomes import annotate_genomes
from kmer_functional_profiler.species import species_index_from_genomes


def test_zero_hit_closed_form() -> None:
    # The plan's worked example: q = 0.95, P_G = 1.
    got = zero_hit_presence(np.array([0.95, 0.95, 1 - 1e-6]), np.array([0.5, 5.0, 50.0]))
    assert np.allclose(got[:2], [0.92, 0.11], atol=0.005)
    assert got[2] < 1e-12  # (almost) always carried, 50 hits expected, none seen: absent


def fake_index(tmp: Path) -> Path:
    """A genome index with only a meta.json: another index for the checksum test."""
    tmp.mkdir(parents=True, exist_ok=True)
    (tmp / "meta.json").write_text(json.dumps({"units": 1000}))
    return tmp


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
    species_index_from_genomes(index_dir, tmp_path / "gi", tmp_path / "c")
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
    other = fake_index(tmp_path / "other")
    with pytest.raises(ValueError, match="another genome index"):
        update(observed([0]), pl.read_csv(tmp / "genomes.tsv", separator="\t"), other, Carriage(c))
