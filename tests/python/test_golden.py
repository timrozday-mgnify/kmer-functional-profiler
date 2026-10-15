"""Amino-acid profiles match the golden outputs of the pre-phase-14 build (tests/golden)."""

from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from kmer_functional_profiler.cli import app
from kmer_functional_profiler.index import read_fasta

ROOT = Path(__file__).resolve().parents[1]
DATA, GOLDEN = ROOT / "data", ROOT / "golden"
READS = [str(DATA / "reads_1.fastq.gz"), str(DATA / "reads_2.fastq.gz")]
INDEXES = {"dense": "--k 7 --t-base 1 --n-min 0", "sparse": "--k 7 --t-base 0.1 --n-min 4",
           "default": ""}  # fmt: skip
QUERIES = {"": ["--aai", "--all-estimators"], "_edges": ["--frames", "edges:30"]}


@pytest.fixture(scope="module")
def members(tmp_path_factory: pytest.TempPathFactory) -> Path:
    rows = list(read_fasta(DATA / "proteins.faa"))
    path = tmp_path_factory.mktemp("golden") / "members.parquet"
    pl.DataFrame(
        {
            "protein_id": [i for i, _ in rows],
            "cluster_rep": [f"u{n // 2}" for n in range(len(rows))],  # pairs share a unit
            "full_length": [True] * len(rows),
            "sequence": [s for _, s in rows],
        }
    ).write_parquet(path)
    return path


@pytest.mark.parametrize("index", INDEXES)
def test_profiles_match_golden(members: Path, tmp_path: Path, index: str) -> None:
    runner = CliRunner()
    built = runner.invoke(
        app, ["index", str(members), str(tmp_path / "idx"), *INDEXES[index].split()]
    )
    assert built.exit_code == 0, built.output
    for suffix, args in QUERIES.items():
        out = tmp_path / f"{index}{suffix}.tsv"
        ran = runner.invoke(app, ["query", str(tmp_path / "idx"), *READS, "--out", str(out), *args])
        assert ran.exit_code == 0, ran.output
        got = pl.read_csv(out, separator="\t", infer_schema_length=None)
        want = pl.read_csv(GOLDEN / out.name, separator="\t", infer_schema_length=None)
        assert got.equals(want), out.name
