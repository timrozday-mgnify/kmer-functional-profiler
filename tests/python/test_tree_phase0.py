"""workflows/tree-index/phase0.py: invariants of the phase-0 measurements on a fixture."""

import json
import os
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "workflows" / "tree-index" / "phase0.py"
SURVEY = SCRIPT.with_name("survey.py")


def run(*args: str, cwd: Path, script: Path = SCRIPT) -> None:
    env = os.environ | {"RAYON_NUM_THREADS": "2", "POLARS_MAX_THREADS": "2"}
    subprocess.run([sys.executable, str(script), *args], cwd=cwd, check=True, env=env)


def test_measure_invariants(tmp_path: Path) -> None:
    run("fixture", "--families", "15", "--members", "8", "--buckets", "1", cwd=tmp_path)
    run(
        "measure", "--members", "members.bucket0.parquet", "--out", "m",
        "--dense-seqs", "100", "--sensitive-seqs", "100", "--chunk-sizes", "50,1000",
        cwd=tmp_path,
    )  # fmt: skip
    m = json.loads((tmp_path / "m" / "measure.json").read_text())
    for tag in ("floors", "no_floors"):
        p = m[tag]["postings"]
        # Fitch is the minimum on its tree; every distinct k-mer costs at least one gain
        assert p["distinct_kmers"] <= p["fitch"] <= min(p["majority"], p["add_only"])
        assert p["fitch"] <= p["flat"]
        assert m[tag]["collapse"]["leaves"] <= m[tag]["collapse"]["distinct_sets"]
    # one level: the sequence-node scheme costs exactly the spanning tree's weight
    flat = m["no_floors"]["postings"]
    assert flat["sequence_mst"] == flat["spanning_tree_weight"]
    # the staged build can only miss edges of the exact spanning tree
    assert m["floors"]["staged"]["staged_weight"] >= m["floors"]["exact_weight"]
    assert all(c["duplicated_postings"] >= 0 for c in m["floors"]["chunks_stored"])


def test_survey_invariants(tmp_path: Path) -> None:
    import polars as pl

    run("fixture", "--families", "15", "--members", "8", "--buckets", "1", cwd=tmp_path)
    run(
        "--members", "members.bucket0.parquet", "--out", "s", "--taus", "0,0.6",
        "--n-singles", "0,16", "--n-1s", "0,1,4",
        cwd=tmp_path, script=SURVEY,
    )  # fmt: skip
    t = pl.read_csv(tmp_path / "s" / "survey.tsv", separator="\t")
    key = ["leaves", "t_base", "n_min", "tau", "n_single"]
    by = t.sort("n_1").group_by(key, maintain_order=True).agg("postings", "blind_proteins", "n_1")
    for postings, blind, n_1 in by.select("postings", "blind_proteins", "n_1").rows():
        assert postings == sorted(postings)  # a larger leaf minimum only adds postings
        # every protein with a valid k-mer is covered once each leaf holds one
        assert all(b == 0 for b, n in zip(blind, n_1, strict=True) if n >= 1)
    # Fitch on the spanning tree can only improve on the sequence-MST scheme (same tree, free
    # internal states); on any tree, Fitch is the cheapest of the tree's schemes
    trees = json.loads((tmp_path / "s" / "survey.json").read_text())["trees"]
    for tree in trees.values():
        for sc in tree["schemes"].values():
            assert sc["fitch_on_spanning_tree"] <= sc["sequence_mst"]
            # the survey's tree is the spanning tree by default
            assert sc["fitch"] == sc["fitch_on_spanning_tree"]
            assert sc["distinct_kmers"] <= sc["fitch"] <= min(sc["majority"], sc["add_only"])
