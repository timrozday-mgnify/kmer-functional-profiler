"""Continuous lineage model prototype (phase 11, step 12): needs the phylo group."""

import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from typer.testing import CliRunner

pytest.importorskip("numpyro")

from kmer_functional_profiler.cli import app  # noqa: E402
from kmer_functional_profiler.lineage import _prune, fit_factors, fit_lineages  # noqa: E402
from kmer_functional_profiler.species import write_species_index  # noqa: E402

CORE, ACCESSORY, E = 30, 25, 10.0


def strain(rng: np.random.Generator, clade: int) -> np.ndarray:
    """Core units, then clade 0's accessory units, then clade 1's: each clade carries its
    own at 0.9 and the other's at 0.1."""
    own = 0.9 if clade == 0 else 0.1
    return np.r_[np.ones(CORE), rng.random(ACCESSORY) < own,
                 rng.random(ACCESSORY) < 1 - own].astype(float)  # fmt: skip


def test_lineage_predicts_accessory_units_from_the_ones_seen_and_splits_two_strains() -> None:
    rng = np.random.default_rng(1)
    carried = np.array([strain(rng, c) for c in [0] * 12 + [1] * 12])
    f = fit_factors(carried)
    e = np.full(carried.shape[1], E)
    q = carried.mean(axis=0)
    # one new strain of clade 0 at 0.05x: most of its units have no hits
    new = strain(rng, 0)
    h = rng.poisson(new * 0.05 * E).astype(float)
    fit = _prune(fit_lineages(h, np.zeros_like(h), e, f), f.sigma)
    assert len(fit["depth"]) == 1
    zero = (h == 0) & (q > 0.1) & (q < 0.9)
    lam = fit["depth"][0]
    species = q * np.exp(-lam * E) / (1 - q + q * np.exp(-lam * E))  # the species model's r

    def loss(p: np.ndarray) -> float:
        p, x = np.clip(p[zero], 1e-6, 1 - 1e-6), new[zero]
        return float(-(x * np.log(p) + (1 - x) * np.log1p(-p)).mean())

    assert loss(fit["carriage"][0]) < 0.9 * loss(species)
    # one strain of each clade at 3x and 2x: two lineages, each nearest its own clade
    other = strain(rng, 1)
    h = rng.poisson((new * 3 + other * 2) * E).astype(float)
    fit = _prune(fit_lineages(h, np.zeros_like(h), e, f), f.sigma)
    assert len(fit["depth"]) == 2
    assert np.allclose(fit["depth"], [3, 2], rtol=0.15)
    nearest = np.linalg.norm(fit["x"][:, None] - f.x[None], axis=-1).argmin(axis=1)
    assert nearest[0] < 12 <= nearest[1]


def test_lineage_cli(tmp_path: Path) -> None:
    rng = np.random.default_rng(2)
    carried = np.array([strain(rng, c) for c in [0] * 8 + [1] * 8])
    (tmp_path / "idx").mkdir()
    (tmp_path / "idx" / "meta.json").write_text("{}")
    g, u = np.nonzero(carried)
    rows = pl.DataFrame({"genome": g, "unit": u, "c": E}).cast(
        {"genome": pl.UInt32, "unit": pl.UInt32}
    )
    genomes = pl.DataFrame(
        {
            "genome": range(len(carried)),
            "name": [f"g{i}" for i in range(len(carried))],
            "taxonomy": "d__B;f__F;g__G;s__S",
        }  # fmt: skip
    ).cast({"genome": pl.UInt32})
    write_species_index(tmp_path / "si", tmp_path / "idx", carried.shape[1], rows, genomes,
                        frozenset(), {})  # fmt: skip
    hits = rng.poisson(strain(rng, 1) * 2 * E)
    pl.DataFrame({"unit": np.arange(len(hits)), "hits": hits}).filter(pl.col("hits") > 0).write_csv(
        tmp_path / "profile.tsv", separator="\t"
    )
    out = [str(tmp_path / f) for f in ("profile.tsv", "si", "lineages.tsv")]
    result = CliRunner().invoke(
        app,
        [
            "lineage",
            *out,
            "--units",
            str(tmp_path / "units.tsv"),
            "--summary",
            str(tmp_path / "summary.json"),
        ],  # fmt: skip
    )
    assert result.exit_code == 0, result.output
    lineages = pl.read_csv(tmp_path / "lineages.tsv", separator="\t")
    assert lineages.height == 1 and lineages["id"][0] == "s__S"
    assert abs(lineages["depth"][0] - 2) < 0.3
    assert int(lineages["nearest"][0].split(",")[0][1:]) >= 8  # a clade-1 genome
    units = pl.read_csv(tmp_path / "units.tsv", separator="\t")
    assert set(units.columns) >= {"species", "lineage", "unit", "carriage_prob", "prevalence"}
    assert json.loads((tmp_path / "summary.json").read_text())["lineages"] == 1
