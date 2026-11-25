"""hmm-annotate's steps on six mini HMMs. The search needs pyhmmer (it runs in kfpfam's image)."""

import gzip
import importlib.util
import math
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import numpy as np
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


def test_seed_settings_for_each_kind_of_threshold(tmp_path: Path) -> None:
    (tmp_path / "models.tsv").write_text(
        "acc\tname\tlength\n"
        + "".join(f"M{i}\tM{i}\t{n}\n" for i, n in enumerate([20, 300, 90, 90]))
    )
    (tmp_path / "t.tsv").write_text(
        "model\tname\tseq\tdom\nM0\tM0\t25.0\t8.0\nM1\tM1\t-inf\t400.0\n"
        "M2\tM2\t21.0\t-inf\nM3\tM3\tinf\tinf\n"
    )
    run(tmp_path, "seeds.py", "bars", "--models", "models.tsv", "--thresholds", "t.tsv",
        "--ga", "ga.tsv", "--gapped", "g.txt")  # fmt: skip
    # the domain's threshold, or the sequence's if that is the only one; capped, and never infinite
    ga = [x.split("\t") for x in (tmp_path / "ga.tsv").read_text().splitlines()]
    assert ga == [[f"M{i}", x, x] for i, x in enumerate(["8.0", "30.0", "21.0", "30.0"])]
    g = list(map(float, (tmp_path / "g.txt").read_text().split()))
    assert g[0] == pytest.approx(
        8 + math.log2(210), abs=0.01
    )  # under the one bar, as a short repeat's is
    assert all(x > 21.25 for x in g[1:])
    run(tmp_path, "seeds.py", "bars", "--models", "models.tsv", "--thresholds", "t.tsv",
        "--ga", "ga.tsv", "--gapped", "g.txt", "--slack", 10)  # fmt: skip
    slack = list(map(float, (tmp_path / "g.txt").read_text().split()))
    assert slack[0] == g[0] and slack[1] > 400 and 21.25 <= slack[2] < g[2]
    # quality falls from 24 to 0 over the grid of T: 9 / .936 is reached at 14.4 steps
    quality = np.arange(24, -1, -1, dtype="<f4")
    np.concatenate([np.zeros(25, "<f4"), quality]).tofile(tmp_path / "s.f32")
    run(tmp_path, "seeds.py", "choose", "--table", "s.f32", "--shape", "10110001", "--theta", 9,
        "--out", "T.txt")  # fmt: skip
    assert float((tmp_path / "T.txt").read_text()) == pytest.approx(
        4.5 + 0.25 * (24 - 9 / 0.936), abs=1e-3
    )


def test_check_summary_totals_and_puts_the_worst_model_first(tmp_path: Path) -> None:
    head = (
        "model\tpairs_truth\tpairs_found\tpairs_extra\tdomains_truth\tdomains_found\tdomains_extra"
    )
    (tmp_path / "c").mkdir()
    (tmp_path / "c/a.tsv").write_text(f"{head}\nA\t4\t4\t0\t4\t4\t0\nB\t1\t1\t0\t2\t1\t0\n")
    (tmp_path / "c/b.tsv").write_text(f"{head}\nB\t2\t1\t0\t2\t1\t0\n")
    run(tmp_path, "tables.py", "check-summary", "c")
    t = pl.read_csv(tmp_path / "check.tsv", separator="\t")
    assert t["model"].to_list() == ["all", "B", "A"]
    assert t["domains_found"].to_list() == [6, 2, 4]
    assert t["domain_recall"].to_list() == [0.75, 0.5, 1.0]


def test_compare_counts_shared_missing_and_extra_domains(tmp_path: Path) -> None:
    here = [
        (10, "PF00001", 1, 50),
        (10, "PF00001", 60, 90),
        (20, "PF00002", 1, 9),
        (30, "PF00001", 1, 9),
    ]
    pl.DataFrame(
        [(*r, True) for r in here] + [(20, "PF00001", 1, 9, False)],
        orient="row", schema=["protein_id", "model", "env_from", "env_to", "passes"],
    ).write_parquet(tmp_path / "d.parquet")  # fmt: skip
    there = [(10, 1, 1, 50), (10, 1, 61, 90), (20, 3, 5, 9), (21, 1, 1, 9), (90, 1, 1, 9)]
    pl.DataFrame(
        there, orient="row", schema=["protein_id", "pfam_accession", "env_from", "env_to"]
    ).write_parquet(tmp_path / "t.parquet")
    run(tmp_path, "compare.py", "--domains", "d.parquet", "--table", "t.parquet", "--lo", 0,
        "--hi", 40, "--sample", 10, "--out", "c.tsv")  # fmt: skip
    # protein 21 is not in the sample and 90 not in the range; PF00001 of 20 does not pass
    assert pl.read_csv(tmp_path / "c.tsv", separator="\t").rows() == [
        ("PF00001", 1, 1, 1, 2, 1, 2), ("PF00002", 0, 0, 1, 0, 0, 1), ("PF00003", 1, 0, 0, 1, 0, 0),
    ]  # fmt: skip


@pytest.mark.parametrize(
    "options", [("--check-sample", 5), ("--all",), ("--all", "--min-cluster", 3)]
)
def test_members_pairs_are_their_representatives(tmp_path: Path, options: tuple) -> None:
    models_tsv(tmp_path)
    thresholds(tmp_path)
    release = PROTEINS.parent
    clusters = pl.read_parquet(release / "mgy_cluster_seqs.parquet")
    members = clusters.filter(pl.col("cluster_member") != pl.col("cluster_rep"))
    reps = clusters["cluster_rep"].unique().sort()
    near = pl.DataFrame(
        {"protein_id": reps, "model": [MODELS[i % 6][0] for i in range(len(reps))]}
    ).filter(pl.col("protein_id") != reps[0])  # one with no hit
    if "--min-cluster" in options:  # the others were not searched, and are members
        small = clusters.filter(pl.len().over("cluster_rep") < 3)["cluster_rep"].unique()
        near = near.filter(~pl.col("protein_id").is_in(small.implode()))
        members = clusters.filter(
            (pl.col("cluster_member") != pl.col("cluster_rep"))
            | pl.col("cluster_rep").is_in(small.implode())
        )
        assert len(small) and near.height
    near.write_parquet(tmp_path / "near.parquet")
    width, check = 40, 5
    cut, end = reps[len(reps) // 2], reps[-1] + 1
    for b, (lo, hi) in enumerate([(0, cut), (cut, end)]):
        run(tmp_path, "members.py", "pairs", "--release", release, "--near", "near.parquet",
            "--lo", lo, "--hi", hi, *options, "--width", width,
            "--out", f"pairs{b}", "--shards", f"shards{b}.txt")  # fmt: skip
    shards = sorted(
        {int(x) for b in "01" for x in (tmp_path / f"shards{b}.txt").read_text().split()}
    )
    seqs = dict(pl.read_parquet(PROTEINS).select("protein_id", "sequence").iter_rows())
    got, names = [], []
    for shard in shards:
        run(tmp_path, "members.py", "fasta", "--release", release, "--pairs", "pairs0", "pairs1",
            "--shard", shard, "--width", width, "--thresholds", "thresholds.tsv",
            "--block", 7, "--out", "m.faa", "--pairs-out", "m.pairs")  # fmt: skip
        lines = (tmp_path / "m.faa").read_text().splitlines()
        ids = [int(x[1:]) for x in lines[::2]]
        assert ids == sorted(ids) and {i // width for i in ids} == {shard}
        assert lines[1::2] == [seqs[i] for i in ids]
        raw = (tmp_path / "m.pairs").read_bytes()
        for at in range(0, len(ids), 7):
            n, residues, size = struct.unpack("<3Q", raw[:24])
            block = lines[1::2][at : at + 7]
            assert (n, residues) == (len(block), sum(map(len, block)))
            pairs = np.frombuffer(raw[24 : 24 + size], "<u4").reshape(-1, 2)
            assert (np.diff(pairs[:, 1].astype(int)) >= 0).all()
            got += [(ids[at + i], MODELS[m][0]) for i, m in pairs]
            raw = raw[24 + size :]
        assert not raw
        names += ids
    candidates = members.join(near, left_on="cluster_rep", right_on="protein_id")
    assert sorted(got) == sorted(candidates.select("cluster_member", "model").rows())
    checked = members.filter((pl.col("cluster_member") % check == 0) | ("--all" in options))
    checked = checked["cluster_member"]
    assert sorted(names) == sorted({*candidates["cluster_member"], *checked})


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

    # near: what passes, and what would with the slack
    t = (tmp_path / "thresholds.tsv").read_text()
    (tmp_path / "thresholds.tsv").write_text(t.replace("Mini1\t25.0\t25.0", "Mini1\t999\t999"))
    for slack in (0, 974):
        searched(tmp_path, seqs, range(6), "--near", "n.parquet", "--near-slack", slack)
        kept = pl.col("model") != "PF90001" if slack == 0 else pl.lit(True)
        near = pl.read_parquet(tmp_path / "n.parquet").rows()
        assert set(near) >= set(labels.filter(kept).rows())
        assert any(m == "PF90001" for _, m in near) == (slack > 0)
    (tmp_path / "thresholds.tsv").write_text(t)

    # every model on every protein, with no prefilter, is the same
    every = searched(tmp_path, seqs, range(0), "--all-pairs")
    assert every.sort(every.columns).equals(d.sort(d.columns))
    stats = pl.read_csv(tmp_path / "s.tsv", separator="\t")
    assert stats.select("proteins", "too_long", "pairs").row(0) == (len(seqs), 0, 6 * len(seqs))

    # so are some models' pairs on stdin and some from a file, with pairs in both
    (tmp_path / "more").write_bytes(stream([len(s) for s in seqs.values()], 30, range(2, 6)))
    both = searched(tmp_path, seqs, range(4), "--more-pairs", "more")
    assert both.equals(d)
    stats = pl.read_csv(tmp_path / "s.tsv", separator="\t")
    assert stats["pairs"].item() == 6 * len(seqs)

    # and nail's seeds, of any order and any number of files, one with none; and those with more
    rows = [f"Mini{m}\t{i}\n" for m in range(1, 7) for i in reversed(seqs)]
    parts = {"a.tsv": rows[:100], "b.tsv.gz": rows[100:], "c.tsv": []}
    for name, part in parts.items():
        data = ("model\tprotein_id\n" + "".join(part)).encode()
        (tmp_path / name).write_bytes(gzip.compress(data) if name.endswith("gz") else data)
    assert searched(tmp_path, seqs, range(0), "--seeds", *parts).equals(d)
    lean = searched(tmp_path, seqs, range(0), "--seeds", "a.tsv", "--block", 30)
    assert 0 < lean.height < d.height
    more = ("--seeds", "a.tsv", "--block", 30, "--more-pairs", "more")
    assert (
        searched(tmp_path, seqs, range(0), *more)
        .filter(pl.col("model") > "PF90002")
        .equals(d.filter(pl.col("model") > "PF90002"))
    )
    for bad in ("Mini1\t999999\n", "Other\t1\n"):
        (tmp_path / "bad.tsv").write_text("model\tprotein_id\n" + bad)
        with pytest.raises(subprocess.CalledProcessError):
            searched(tmp_path, seqs, range(0), "--seeds", "bad.tsv")
    d = searched(tmp_path, seqs, range(6))

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

    # two sets of pairs in different blocks are an error
    (tmp_path / "more").write_bytes(stream([len(s) for s in seqs.values()], 20, range(6)))
    with pytest.raises(subprocess.CalledProcessError):
        searched(tmp_path, seqs, range(4), "--more-pairs", "more")


@needs_pyhmmer
def test_proteins_too_long_for_hmmer_are_counted_not_searched(tmp_path: Path) -> None:
    models_tsv(tmp_path)
    thresholds(tmp_path)
    run(tmp_path, "annotate.py", "press", "--hmm", HMM, "--models", "models.tsv", "--out", "m")
    seqs = dict(pl.read_parquet(PROTEINS).select("protein_id", "sequence").iter_rows())
    hit = searched(tmp_path, seqs, range(6))["protein_id"][0]
    filler = "GSTNQDEKRA"
    long = {1: filler * 6_800 + seqs[hit], 2: filler * 10_000 + seqs[hit], 3: "M" * 50}
    assert set(searched(tmp_path, long, range(6))["protein_id"]) == {1}
    assert pl.read_csv(tmp_path / "s.tsv", separator="\t")["too_long"].item() == 1
    assert set(searched(tmp_path, long, range(0), "--all-pairs")["protein_id"]) == {1}
    assert pl.read_csv(tmp_path / "s.tsv", separator="\t")["too_long"].item() == 1
