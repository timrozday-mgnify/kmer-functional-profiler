"""hmm-annotate's steps on six mini HMMs. The search needs pyhmmer (it runs in kfpfam's image)."""

import importlib.util
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "workflows" / "hmm-annotate"
HMM = ROOT / "tests" / "data" / "mini_fmh" / "Pfam-mini.hmm"
PROTEINS = ROOT / "tests" / "data" / "mini_reannotate" / "mgy_protein_sequences.parquet"
needs_pyhmmer = pytest.mark.skipif(
    importlib.util.find_spec("pyhmmer") is None or shutil.which("hmmsearch") is None,
    reason="needs pyhmmer and HMMER",
)
MODELS = [(f"PF{90000 + i}", f"Mini{i}") for i in range(1, 7)]


def run(cwd: Path, script: str, *args: object, stdin: bytes = b"") -> None:
    subprocess.run([sys.executable, WF / script, *map(str, args)], cwd=cwd, check=True, input=stdin)


def models_tsv(cwd: Path, models: list[tuple[str, str]] = MODELS) -> None:
    (cwd / "models.tsv").write_text("acc\tname\n" + "".join(f"{a}\t{n}\n" for a, n in models))


def thresholds(cwd: Path, *args: object) -> list[tuple[str, float, float]]:
    run(cwd, "tables.py", "thresholds", "--hmm", HMM, "--models", "models.tsv", *args)
    t = pl.read_csv(cwd / "thresholds.tsv", separator="\t")
    return list(t.select("model", "seq", "dom").iter_rows())


def test_thresholds_from_the_hmm_file_and_from_a_table(tmp_path: Path) -> None:
    models_tsv(tmp_path)
    assert thresholds(tmp_path) == [(a, 25.0, 25.0) for a, _ in MODELS]
    # KOfam's ko_list: further columns, names for models with no accession, '-' for no threshold
    models_tsv(tmp_path, [(n, n) for _, n in MODELS])
    rows = ["knum\tthreshold\tscore_type\tprofile_type"]
    rows += [f"Mini{i}\t{10 * i}\t{'full' if i % 2 else 'domain'}\tall" for i in range(1, 6)]
    (tmp_path / "ko_list").write_text("\n".join([*rows, "Mini6\t-\t-\t-"]) + "\n")
    inf = float("inf")
    assert thresholds(tmp_path, "--table", "ko_list") == [
        ("Mini1", 10, -inf), ("Mini2", -inf, 20), ("Mini3", 30, -inf), ("Mini4", -inf, 40),
        ("Mini5", 50, -inf), ("Mini6", inf, inf),
    ]  # fmt: skip


def test_thresholds_fail_on_a_model_with_none(tmp_path: Path) -> None:
    models_tsv(tmp_path)
    lines = HMM.read_text().splitlines(keepends=True)
    (tmp_path / "no_ga.hmm").write_text("".join(x for x in lines if not x.startswith("GA ")))
    with pytest.raises(subprocess.CalledProcessError):
        run(tmp_path, "tables.py", "thresholds", "--hmm", "no_ga.hmm", "--models", "models.tsv")
    (tmp_path / "t.tsv").write_text("model\tthreshold\tscore_type\nPF90001\t1\tfull\n")
    with pytest.raises(subprocess.CalledProcessError):
        thresholds(tmp_path, "--table", "t.tsv")


def test_check_summary_totals_and_puts_the_worst_model_first(tmp_path: Path) -> None:
    head = (
        "model\tpairs_truth\tpairs_found\tpairs_extra\tdomains_truth\tdomains_found\tdomains_extra"
    )
    (tmp_path / "a.tsv").write_text(f"{head}\nA\t4\t4\t0\t4\t4\t0\nB\t1\t1\t0\t2\t1\t0\n")
    (tmp_path / "b.tsv").write_text(f"{head}\nB\t2\t1\t0\t2\t1\t0\n")
    run(tmp_path, "tables.py", "check-summary", "a.tsv", "b.tsv")
    t = pl.read_csv(tmp_path / "check.tsv", separator="\t")
    assert t["model"].to_list() == ["all", "B", "A"]
    assert t["domains_found"].to_list() == [6, 2, 4]
    assert t["domain_recall"].to_list() == [0.75, 0.5, 1.0]


def stream(lengths: list[int], block: int, models: range) -> bytes:
    """Every (protein, model) pair, as `filter join` writes them."""
    out = b""
    for at in range(0, len(lengths), block):
        n = len(lengths[at : at + block])
        pairs = [(i, m) for m in models for i in range(n)]
        data = b"".join(struct.pack("<II", *p) for p in pairs)
        out += struct.pack("<3Q", n, sum(lengths[at : at + block]), len(data)) + data
    return out


def searched(cwd: Path, seqs: dict[int, str], models: range, *args: object) -> pl.DataFrame:
    (cwd / "p.faa").write_text("".join(f">{i}\n{s}\n" for i, s in seqs.items()))
    run(
        cwd, "annotate.py", "search", "--faa", "p.faa", "--pressed", "m", "--threads", 2,
        "--thresholds", "thresholds.tsv", "--domains", "d.parquet", "--labels", "l.parquet",
        "--stats", "s.tsv", *args, stdin=stream([len(s) for s in seqs.values()], 30, models),
    )  # fmt: skip
    return pl.read_parquet(cwd / "d.parquet")


@needs_pyhmmer
def test_search_gives_hmmsearch_cut_ga(tmp_path: Path) -> None:
    models_tsv(tmp_path)
    thresholds(tmp_path)
    run(tmp_path, "annotate.py", "press", "--hmm", HMM, "--models", "models.tsv", "--out", "m")
    seqs = dict(pl.read_parquet(PROTEINS).select("protein_id", "sequence").iter_rows())
    d = searched(tmp_path, seqs, range(6))
    subprocess.run(
        ["hmmsearch", "--cut_ga", "--domtblout", "ga.tbl", "-o", "/dev/null", HMM, "p.faa"],
        cwd=tmp_path,
        check=True,
    )
    expected = sorted(
        (int(c[0]), c[4].split(".")[0], *map(int, c[15:21]), round(float(c[13]), 1))
        for c in map(str.split, (tmp_path / "ga.tbl").read_text().splitlines())
        if not c[0].startswith("#")
    )
    cols = ["protein_id", "model", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to"]
    got = sorted((*r[:-1], round(r[-1], 1)) for r in d.select(*cols, "score").iter_rows())
    assert expected and [g[:-1] for g in got] == [e[:-1] for e in expected]
    assert all(abs(g[-1] - e[-1]) < 0.11 for g, e in zip(got, expected, strict=True))
    assert d["passes"].all()
    for r in d.iter_rows(named=True):
        seq = seqs[r["protein_id"]]
        assert r["left"] == min(10, r["ali_from"] - 1)
        assert r["right"] == min(10, len(seq) - r["ali_to"])
        assert r["residues"] == seq[r["ali_from"] - 1 - r["left"] : r["ali_to"] + r["right"]]
    labels = pl.read_parquet(tmp_path / "l.parquet")
    assert sorted(labels.rows()) == sorted(set(d.select("protein_id", "model").rows()))

    # the check finds it whole, and finds what is missing when a model is left out
    check = (
        "annotate.py", "check", "--faa", "p.faa", "--hmm", HMM, "--thresholds", "thresholds.tsv",
        "--domains", "d.parquet", "--out", "c.tsv",
    )  # fmt: skip
    run(tmp_path, *check)
    c = pl.read_csv(tmp_path / "c.tsv", separator="\t")
    assert c["domains_truth"].sum() == c["domains_found"].sum() == d.height
    assert c["domains_extra"].sum() == 0
    searched(tmp_path, seqs, range(1, 6))
    run(tmp_path, *check)
    c = pl.read_csv(tmp_path / "c.tsv", separator="\t").filter(pl.col("model") == "PF90001")
    assert c["domains_truth"].item() > 0 and c["domains_found"].item() == 0

    # a model with no threshold never passes, but shows with --weak-evalue
    t = (tmp_path / "thresholds.tsv").read_text()
    (tmp_path / "thresholds.tsv").write_text(t.replace("Mini1\t25.0\t25.0", "Mini1\tinf\tinf"))
    one = pl.col("model") == "PF90001"
    assert searched(tmp_path, seqs, range(6)).filter(one).is_empty()
    weak = searched(tmp_path, seqs, range(6), "--weak-evalue", 1e-3).filter(one)
    assert weak.height >= 3 and not weak["passes"].any()


@needs_pyhmmer
def test_proteins_too_long_for_the_filter_or_hmmer(tmp_path: Path) -> None:
    models_tsv(tmp_path)
    thresholds(tmp_path)
    run(tmp_path, "annotate.py", "press", "--hmm", HMM, "--models", "models.tsv", "--out", "m")
    seqs = dict(pl.read_parquet(PROTEINS).select("protein_id", "sequence").iter_rows())
    d = searched(tmp_path, seqs, range(6))
    hit = d["protein_id"][0]
    filler = "GSTNQDEKRA"
    long = {1: filler * 6_800 + seqs[hit], 2: filler * 10_000 + seqs[hit], 3: "M" * 50}
    # no pairs from the filter: the long protein is found all the same, the longest is counted
    d = searched(tmp_path, long, range(0))
    assert set(d["protein_id"]) == {1}
    assert pl.read_csv(tmp_path / "s.tsv", separator="\t")["too_long"].item() == 1
