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
         "--frames", "stopfree", "edges", "--out", str(tmp_path)],
        check=True,
        capture_output=True,
    )  # fmt: skip
    summary = pl.read_csv(tmp_path / "summary.tsv", separator="\t")
    assert summary.height == 16  # one row per rule and frame mode
    assert (summary["completeness"] > 0.5).all()
    zi = summary.filter(rule="gather_zi")
    assert zi["aai_n"].min() > 0  # type: ignore[operator]
    assert zi["aai_bias_1.0"].abs().max() < 0.05  # type: ignore[operator]
    assert {"complete_200", "bias_len_450"} <= set(summary.columns)
