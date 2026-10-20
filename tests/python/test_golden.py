"""Amino-acid profiles match the golden outputs of the pre-phase-14 build (tests/golden): exact
for counts and ids, to a relative 1e-6 for floats (EM and fits differ in the last digits
between platforms; the plan's rule for golden tests)."""

from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal
from typer.testing import CliRunner

from kmer_functional_profiler import index as index_module
from kmer_functional_profiler import query
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


@pytest.mark.parametrize("in_order", [False, True])
@pytest.mark.parametrize("index", INDEXES)
def test_profiles_match_golden(
    members: Path, tmp_path: Path, index: str, in_order: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    if in_order:  # tier 1 and unit columns read in small blocks, never in place (step 40)
        monkeypatch.setattr(index_module, "PAGE", 2**62)  # every block read whole
        monkeypatch.setattr(index_module, "READ_ROWS", 64)
        monkeypatch.setattr(query, "WALK_BUCKETS", 4)
        monkeypatch.setattr(query, "WALK_HASHES", 16)
        monkeypatch.setattr(query, "AAI_CHUNK", 3)
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
        # present_q (2026-10-20) is new beside the frozen outputs: the q-value of present_prob
        q = query.fdr_q(got["present_prob"].to_numpy())
        assert got["present_q"].to_numpy() == pytest.approx(q)
        got = got.drop("present_q")
        want = pl.read_csv(GOLDEN / out.name, separator="\t", infer_schema_length=None)
        if "--aai" not in args:  # no AAI: the detected units only, labelled within them
            want = (
                want.filter(pl.col("gather_rank").is_not_null())
                .drop("aai_naive", "aai_naive_lower_bound")
                .with_columns(component=pl.col("unit").min().over("component"))
            )
        assert_frame_equal(got, want, check_exact=False, rel_tol=1e-6, abs_tol=1e-9)
