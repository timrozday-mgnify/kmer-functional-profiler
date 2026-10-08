"""The study-index recipe's tables, and `index --like`."""

import importlib.util
import json
from pathlib import Path

from typer.testing import CliRunner

from kmer_functional_profiler.cli import app

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("study", ROOT / "workflows/study-index/study.py")
assert SPEC is not None and SPEC.loader is not None
study = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(study)
SEQS = ["MKVLAAGIVGLLLAAPAQA", "MKVLAAGIVGLLLAAPAQS", "MTEYKLVVVGAGGVGKSAL"]


def test_members_and_index_like(tmp_path: Path) -> None:
    faa = tmp_path / "genes.faa"
    faa.write_text(
        f">a_1 # 1 # 57 # 1 # ID=1_1;partial=00\n{SEQS[0]}*\n"
        f">a_2 # 1 # 57 # 1 # ID=1_2;partial=10\n{SEQS[1]}\n"
        f">b_1\n{SEQS[2]}\n"
    )
    clusters = tmp_path / "clusters.tsv"
    clusters.write_text("a_1\ta_1\na_1\ta_2\n")  # b_1 left out: becomes a singleton
    members = study.members([str(faa)], str(clusters))
    assert members.rows() == [
        ("a_1", "a_1", True, SEQS[0]),
        ("a_2", "a_1", False, SEQS[1]),
        ("b_1", "b_1", True, SEQS[2]),
    ]
    domtbl = tmp_path / "x.domtbl"
    domtbl.write_text("# header\n" + " ".join(["a_1", "-", "1", "PF1", "PF00001.20"] + ["0"] * 18))
    assert study.pfam([str(domtbl)]).rows() == [("a_1", "PF00001")]

    members.write_parquet(tmp_path / "m.parquet")
    runner = CliRunner()
    base = ["index", str(tmp_path / "m.parquet"), str(tmp_path / "base")]
    assert (
        runner.invoke(app, [*base, "--k", "5", "--t-base", "0.5", "--t-dense", "1"]).exit_code == 0
    )
    like = ["index", str(tmp_path / "m.parquet"), str(tmp_path / "like"), "--k", "7"]
    assert runner.invoke(app, [*like, "--like", str(tmp_path / "base")]).exit_code == 0
    params = [
        json.loads((tmp_path / d / "meta.json").read_text())["params"] for d in ("base", "like")
    ]
    assert params[0] == params[1] and params[1]["k"] == 5
