"""workflows/presence-calibration/calibrate.py: score and decide on hand-made profiles."""

import argparse
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflows" / "presence-calibration"))
import calibrate  # noqa: E402


def _profile(
    path: Path,
    q: dict[int, float],
    component: dict[int, int] | None = None,
    kmers: dict[int, tuple[int, int]] | None = None,  # unit -> (kmers_hit, kmers_unique)
) -> Path:
    comp, km = component or {}, kmers or {}
    pl.DataFrame(
        {
            "unit": list(q),
            "present_q": list(q.values()),
            "component": [comp.get(u, u) for u in q],
            "kmers_hit": [km.get(u, (3, 3))[0] for u in q],
            "kmers_unique": [km.get(u, (3, 3))[1] for u in q],
        }
    ).write_csv(path, separator="\t")
    return path


def test_score_counts_spike_recall_new_false_and_null(tmp_path: Path) -> None:
    pl.DataFrame(
        {"unit": [1, 2, 3, 4], "identity_bin": [1.0, 1.0, 0.9, 0.9], "depth": [0.1, 3.0, 3.0, 10.0]}
    ).write_csv(tmp_path / "truth.tsv", separator="\t")
    pl.DataFrame({"unit": [1, 2, 3, 4], "uniq": [0.1, 0.9, 0.9, 0.9]}).write_csv(
        tmp_path / "units.tsv", separator="\t"
    )
    # spiked: 2, 3 called at 0.05 (4 only at 0.1); 7 already in the unspiked sample; 8 new,
    # false, in spiked 3's component (a relative); 10 new, false, unrelated. Spiked 1 is hit
    # (one shared k-mer) but not gathered; 4 has two hit k-mers, one its own.
    spiked = _profile(
        tmp_path / "p.tsv",
        {1: 1.0, 2: 0.0, 3: 0.04, 4: 0.08, 7: 0.0, 8: 0.01, 9: 0.5, 10: 0.03},
        component={3: 3, 8: 3},
        kmers={1: (1, 0), 4: (2, 1)},
    )
    unspiked = _profile(tmp_path / "u.tsv", {7: 0.0, 8: 0.5})
    null = _profile(tmp_path / "n.tsv", {5: 0.02})
    out = tmp_path / "scores.tsv"
    calibrate.score(
        argparse.Namespace(
            sample="seed1", share="0.1", truth=tmp_path / "truth.tsv",
            units=tmp_path / "units.tsv", profile=spiked, unspiked=unspiked, null=null, out=out,
        )
    )  # fmt: skip
    s = pl.read_csv(out, separator="\t").filter(pl.col("q") == 0.05).row(0, named=True)
    assert (s["called"], s["new_false"], s["unspiked_called"], s["null_called"]) == (5, 2, 1, 1)
    assert s["recall"] == 0.5
    assert s["false_share"] == 2 / 4
    assert s["new_relative"] == 1 and s["false_share_unrelated"] == 1 / 3
    assert (s["recall_hit1"], s["recall_hit2"], s["recall_gathered"]) == (1.0, 0.75, 0.75)
    assert s["recall_hit1_depth_0"] == 1.0 and s["recall_hit1_identity_0.9"] == 1.0
    assert (s["unspiked_hit1"], s["unspiked_hit2"], s["null_hit1"]) == (2, 2, 1)
    assert s["recall_depth_0"] == 0.0 and s["recall_depth_2"] == 1.0 and s["recall_depth_5"] == 0.0
    assert s["recall_identity_0.9"] == 0.5 and s["n_uniq_0.75"] == 3

    # decide: share 0.1 fails on false_share; 0.01 passes, so it is chosen
    other = pl.read_csv(out, separator="\t").with_columns(
        share=pl.lit(0.01), new_false=pl.lit(0), false_share=pl.lit(0.0), null_called=pl.lit(0)
    )
    other.write_csv(tmp_path / "s2.tsv", separator="\t")
    calibrate.decide(
        argparse.Namespace(
            scores=[out, tmp_path / "s2.tsv"], q=0.05, max_false=0.05, max_null=0.01,
            out=tmp_path / "summary.tsv", decision=tmp_path / "decision.tsv",
        )
    )  # fmt: skip
    d = pl.read_csv(tmp_path / "decision.tsv", separator="\t")
    assert d.filter(pl.col("chosen"))["share"].to_list() == [0.01]
    assert (tmp_path / "scores.tsv").exists()
    summary = pl.read_csv(tmp_path / "summary.tsv", separator="\t")
    assert {"recall_hit1", "unspiked_hit2", "null_hit1"} <= set(summary.columns)


def test_reverse_reverses_sequence_and_quality(tmp_path: Path) -> None:
    import gzip

    for m in (1, 2):
        with gzip.open(tmp_path / f"r{m}.fq.gz", "wt") as f:
            f.write("@a\nACGT\n+\nIIJK\n")
    calibrate.reverse(
        argparse.Namespace(
            reads=[tmp_path / "r1.fq.gz", tmp_path / "r2.fq.gz"],
            out=[tmp_path / "n1.fq.gz", tmp_path / "n2.fq.gz"],
        )
    )
    with gzip.open(tmp_path / "n1.fq.gz", "rt") as f:
        assert f.read() == "@a\nTGCA\n+\nKJII\n"


def test_back_translate_takes_ambiguous_residues() -> None:
    import random

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflows" / "spike-benchmark"))
    import spike

    cds = spike.back_translate("MXBZJUOK", random.Random(0))
    protein = spike.translate(cds)  # final stop dropped, internal stops as X
    assert len(cds) == 27 and cds[-3:] in spike.SYNONYMS["*"]
    assert protein[0] == "M" and protein[-1] == "K" and "X" not in protein
