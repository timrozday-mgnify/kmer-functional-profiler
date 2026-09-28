"""The simulated abundance benchmark runs end to end on a tiny reference."""

import subprocess
import sys
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflows" / "sim-benchmark" / "sim.py"


def test_sim_benchmark(tmp_path: Path) -> None:
    subprocess.run(
        [sys.executable, SCRIPT, "--families", "10", "--seeds", "1", "--configs", "dense",
         "--seed-fasta", str(ROOT / "tests" / "data" / "proteins.faa"),
         "--out", str(tmp_path)],
        check=True,
        capture_output=True,
    )  # fmt: skip
    summary = pl.read_csv(tmp_path / "summary.tsv", separator="\t")
    assert summary.height == 7  # one row per rule
    assert (summary["completeness"] > 0.5).all()
