"""bundle.py on a made-up family: the right label comes back, and what should be refused is. Run it:
test_bundle.py"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import polars as pl

BUNDLE = Path(__file__).resolve().parent / "bundle.py"
AA = "ACDEFGHIKLMNPQRSTVWY"


def run(*a, ok=True, says="") -> str:
    r = subprocess.run([sys.executable, BUNDLE, *map(str, a)], capture_output=True, text=True)
    assert (r.returncode == 0) == ok, (a, r.returncode, r.stderr[-800:])
    assert says in r.stdout + r.stderr, (says, r.stderr[-800:])
    return r.stdout + r.stderr


with tempfile.TemporaryDirectory() as tmp:
    d = Path(tmp)
    rng = np.random.default_rng(0)
    rand = lambda n: "".join(rng.choice(list(AA), n))  # noqa: E731
    (d / "m.hmm").write_text("ACC   PF00007.1\nLENG  80\n//\nACC   PF00009.3\nLENG  30\n//\n")
    (d / "other.hmm").write_text("ACC   PF00007.2\nLENG  90\n//\n")
    core = rand(80)
    dom = {
        i: "".join(c if rng.random() < 0.6 else rng.choice(list(AA)) for c in core)
        for i in range(1, 31)
    }
    pre, post = {i: rand(i % 14) for i in dom}, {i: rand(i % 5) for i in dom}
    # protein 30 has an insert of 4 after column 20, and no column 41-45
    ins = rand(4)
    body = {i: s for i, s in dom.items()}
    body[30] = dom[30][:20] + ins + dom[30][20:40] + dom[30][45:]
    cig = {i: "80M" for i in dom} | {30: "20M4I20M5D35M"}
    pl.DataFrame(
        [(i, 7, 1, len(pre[i]) + 1, len(pre[i]) + len(body[i]), cig[i], i != 29) for i in dom],
        orient="row",
        schema=[
            "protein_id",
            "pfam_accession",
            "hmm_from",
            "ali_from",
            "ali_to",
            "cigar",
            "above_ga",
        ],
    ).write_parquet(d / "h.parquet")
    pl.DataFrame(
        [(i + 5 * 10**9 * (i == 30), pre[i] + body[i] + post[i]) for i in dom],
        orient="row",
        schema=["protein_id", "sequence"],
    ).write_parquet(d / "s.parquet")
    # protein 30's accession needs 5 bytes
    h = pl.read_parquet(d / "h.parquet").with_columns(
        protein_id=pl.col("protein_id") + 5 * 10**9 * (pl.col("protein_id") == 30)
    )
    h.write_parquet(d / "h.parquet")
    (d / "ko.tsv").write_text(
        "".join(f"{i}\tK{i:05d}\n" for i in dom if i % 2 == 0 and i < 30)
        + "5000000030\tK00030\n5000000030\tK99999\n"
    )
    args = [
        "--hits",
        d / "h.parquet",
        "--seqs",
        d / "s.parquet",
        "--hmm",
        d / "m.hmm",
        "--label",
        f"ko={d / 'ko.tsv'}",
        "--pfam",
        "t",
    ]
    run("build", *args, "--out", d / "b", "--collapse", "0.7")
    run("verify", d / "b", says="as built")
    m = json.loads((d / "b/manifest.json").read_text())
    assert m["rows"] == 29, m["rows"]  # 29 is under GA; unlabelled rows stay
    reads = [
        ("r2", 7, 5, "40M", dom[2][4:44]),
        ("r3", 7, 5, "40M", dom[3][4:44]),
        ("r30", 7, 11, "10M4I20M5D10M", dom[30][10:20] + ins + dom[30][20:40] + dom[30][45:55]),
        ("r30", 7, 60, "20M", rand(20)),
        ("junk", 9, 1, "30M", rand(30)),
    ]
    s = pl.DataFrame(
        reads, orient="row", schema=["read", "f", "hmm_from", "cigar", "residues"]
    ).with_columns(pfam=pl.format("PF{}", pl.col("f").cast(pl.String).str.zfill(5)))
    s.drop("f").write_parquet(d / "x.parquet")
    out = run(
        "profile",
        "--bundle",
        d / "b",
        "--parquet",
        d / "x.parquet",
        "--out",
        d / "p",
        says="warning",
    )
    c = {
        r["read"]: r
        for r in pl.read_parquet(d / "p.reads.parquet")
        .sort("score", nulls_last=False)
        .iter_rows(named=True)
    }
    assert (
        c["r2"]["ko"] == "K00002"
        and c["r2"]["protein"] == "MGYP000000000002"
        and c["r2"]["similarity"] == 1
    )
    assert c["r3"]["ko"] is None and c["r3"]["protein"] == "MGYP000000000003"
    assert c["r30"]["ko"] == "K00030+K99999" and c["r30"]["protein"] == "MGYP005000000030", c["r30"]
    assert c["junk"]["score"] is None  # a family of the models with no row in the index
    p = dict(pl.read_csv(d / "p.ko.tsv", separator="\t").select("label", "reads").iter_rows())
    assert p == {"K00002": 1, "K00030+K99999": 1}, p  # r30 counts once
    s = json.loads((d / "p.summary.json").read_text())
    assert (s["reads"], s["reads_called"], s["labelled_reads"]["ko"]) == (4, 3, 2), s

    run("build", *args, "--out", d / "b", ok=False, says="not empty")
    run(
        "profile",
        "--bundle",
        d / "b",
        "--parquet",
        d / "x.parquet",
        "--out",
        d / "p",
        "--hmm",
        d / "m.hmm",
    )
    run(
        "profile",
        "--bundle",
        d / "b",
        "--parquet",
        d / "x.parquet",
        "--out",
        d / "p",
        "--hmm",
        d / "other.hmm",
        ok=False,
        says="not the bundle's Pfam",
    )
    # the pipeline's version column: the bundle's models pass, others do not
    x = pl.read_parquet(d / "x.parquet").with_columns(
        above_ga=pl.col("read") != "r3",
        pfam_version=pl.when(pl.col("pfam") == "PF00007").then(1).otherwise(3).cast(pl.Int16),
    )
    x.write_parquet(d / "z.parquet")
    assert "warning" not in run(
        "profile", "--bundle", d / "b", "--parquet", d / "z.parquet", "--out", d / "p"
    )
    assert json.loads((d / "p.summary.json").read_text())["reads"] == 4
    run("profile", "--bundle", d / "b", "--parquet", d / "z.parquet", "--out", d / "p", "--ga-only")
    assert json.loads((d / "p.summary.json").read_text())["reads"] == 3
    x.with_columns(pfam_version=pl.lit(2, pl.Int16)).write_parquet(d / "y.parquet")
    run(
        "profile",
        "--bundle",
        d / "b",
        "--parquet",
        d / "y.parquet",
        "--out",
        d / "p",
        ok=False,
        says="other versions",
    )
    # past the model's end, as the reads of another release can be
    pl.read_parquet(d / "x.parquet").with_columns(hmm_from=pl.lit(60)).write_parquet(
        d / "w.parquet"
    )
    run(
        "profile",
        "--bundle",
        d / "b",
        "--parquet",
        d / "w.parquet",
        "--out",
        d / "p",
        ok=False,
        says="another Pfam release?",
    )
    pl.read_parquet(d / "x.parquet").with_columns(pfam=pl.lit("PF00008")).write_parquet(
        d / "v.parquet"
    )
    run(
        "profile",
        "--bundle",
        d / "b",
        "--parquet",
        d / "v.parquet",
        "--out",
        d / "p",
        ok=False,
        says="not in the bundle's Pfam",
    )
    run(
        "build",
        *[x if x != d / "m.hmm" else d / "other.hmm" for x in args],
        "--out",
        d / "b2",
        "--collapse",
        "0",
    )  # PF00007 alone: still builds
    h.with_columns(pfam_accession=pl.lit(8)).write_parquet(d / "h8.parquet")
    run(
        "build",
        *[x if x != d / "h.parquet" else d / "h8.parquet" for x in args],
        "--out",
        d / "b3",
        ok=False,
        says="another release?",
    )
    # a byte changed, then a file cut
    f = d / "b/rows.bin"
    os.chmod(f, 0o644)
    raw = bytearray(f.read_bytes())
    raw[-3] ^= 1
    f.write_bytes(raw)
    run("verify", d / "b", ok=False, says="changed since it was built: rows.bin")
    f.write_bytes(raw[:-1])
    run(
        "profile",
        "--bundle",
        d / "b",
        "--parquet",
        d / "z.parquet",
        "--out",
        d / "p",
        ok=False,
        says="missing, or not",
    )

    # members: a read of member 102 of cluster 101, whose representative is far from it while
    # cluster 201's is near. Representatives alone call 201; with members, 101. Member 103 is
    # the same row as 102 (kept once); 302's representative has no row here (left out).
    x = rand(80)
    mut = lambda s, p: "".join(c if rng.random() > p else rng.choice(list(AA)) for c in s)  # noqa: E731
    doms = {101: mut(x, 0.5), 102: x, 103: x, 201: mut(x, 0.15), 302: mut(x, 0.3)}
    cluster = {101: 101, 102: 101, 103: 101, 201: 201, 302: 301}
    flank = rand(12)
    pl.DataFrame(
        [(i, 11, 1, 13, 92, "80M", True) for i in doms],
        orient="row",
        schema=[
            "protein_id",
            "pfam_accession",
            "hmm_from",
            "ali_from",
            "ali_to",
            "cigar",
            "above_ga",
        ],
    ).write_parquet(d / "mh.parquet")
    pl.DataFrame(
        [(i, flank + s + flank) for i, s in doms.items()] + [(301, rand(100))],
        orient="row",
        schema=["protein_id", "sequence"],
    ).write_parquet(d / "ms.parquet")
    pl.DataFrame(
        list(cluster.items()) + [(301, 301)], orient="row", schema=["protein_id", "cluster_rep"]
    ).write_parquet(d / "mc.parquet")
    (d / "mko.tsv").write_text("102\tK00102\n103\tK00102\n")
    (d / "mec.tsv").write_text("102\t1.1.1.1\n103\t1.1.1.1\n")
    (d / "mm.hmm").write_text("ACC   PF00011.1\nLENG  80\n//\n")
    margs = [
        "--hits",
        d / "mh.parquet",
        "--seqs",
        d / "ms.parquet",
        "--hmm",
        d / "mm.hmm",
        "--label",
        f"ko={d / 'mko.tsv'}",
        "--label",
        f"ec={d / 'mec.tsv'}",  # not in alphabetical order: calls must not swap the schemes
    ]
    run(
        "build",
        *margs,
        "--clusters",
        d / "mc.parquet",
        "--out",
        d / "bm",
        ok=False,
        says="--collapse 0",
    )
    run("build", *margs, "--clusters", d / "mc.parquet", "--collapse", "0", "--out", d / "bm")
    # representatives only, as an index without members holds them
    reps = [k for k, v in cluster.items() if k == v]
    pl.read_parquet(d / "mh.parquet").filter(pl.col("protein_id").is_in(reps)).write_parquet(
        d / "mr.parquet"
    )
    run(
        "build",
        *[d / "mr.parquet" if x == d / "mh.parquet" else x for x in margs],
        "--collapse",
        "0",
        "--out",
        d / "br",
    )
    m = json.loads((d / "bm/manifest.json").read_text())
    assert (m["rows"], m["member_rows"], m["member_rows_unreached"]) == (2, 1, 1), m
    pl.DataFrame(
        [("m", "PF00011", 6, "50M", x[5:55])],
        orient="row",
        schema=["read", "pfam", "hmm_from", "cigar", "residues"],
    ).write_parquet(d / "mx.parquet")
    for b, want, ko in (("bm", "MGYP000000000101", "K00102"), ("br", "MGYP000000000201", None)):
        run(
            "profile",
            "--bundle",
            d / b,
            "--parquet",
            d / "mx.parquet",
            "--out",
            d / "mp",
            "--hmm",
            d / "mm.hmm",
        )
        r = pl.read_parquet(d / "mp.reads.parquet").row(0, named=True)
        assert (r["protein"], r["ko"], r["ties"]) == (want, ko, 1), (b, r)
        assert r["ec"] == (ko and "1.1.1.1"), (b, r)

    # in parts: the members fixture as two families, each built on its own and then assembled,
    # gives each family's block byte for byte, and the same calls, as one build of both
    sys.path.insert(0, str(BUNDLE.parent))
    import argparse

    import bundle

    h2 = pl.read_parquet(d / "mh.parquet")
    pl.concat(
        [h2, h2.with_columns(pfam_accession=pl.lit(12, h2["pfam_accession"].dtype))]
    ).write_parquet(d / "mh2.parquet")
    (d / "mm2.hmm").write_text("ACC   PF00011.1\nLENG  80\n//\nACC   PF00012.1\nLENG  80\n//\n")
    one = argparse.Namespace(
        hits=d / "mh2.parquet",
        seqs=d / "ms.parquet",
        clusters=d / "mc.parquet",
        label=[f"ko={d / 'mko.tsv'}"],
    )
    rows, _ = bundle.rows_of(one, {11: 80, 12: 80})
    for f in (11, 12):
        rows.filter(pl.col("pfam_accession") == f).write_parquet(d / f"part{f}.parquet")
    build2 = ["--hits", d / "mh2.parquet", "--seqs", d / "ms.parquet", "--hmm", d / "mm2.hmm"]
    build2 += ["--label", f"ko={d / 'mko.tsv'}", "--clusters", d / "mc.parquet", "--collapse", "0"]
    run("build", *build2, "--out", d / "whole")
    for f in (11, 12):
        run(
            "build-part",
            "--rows",
            d / f"part{f}.parquet",
            "--hmm",
            d / "mm2.hmm",
            "--out",
            d / f"p{f}",
        )
    run(
        "assemble",
        "--parts",
        d / "p11",
        d / "p12",
        "--hmm",
        d / "mm2.hmm",
        "--out",
        d / "joined",
    )
    run("verify", d / "joined", says="as built")

    def blocks(b: Path) -> dict[int, bytes]:
        raw = (b / "rows.bin").read_bytes()
        out = {}
        for line in (b / "families.tsv").read_text().splitlines():
            c = [int(x) for x in line.split("\t")]
            out[c[0]] = raw[c[4] : c[4] + c[5]]
        return out

    assert blocks(d / "whole") == blocks(d / "joined") and len(blocks(d / "joined")) == 2
    labels = (d / "whole/labels.ko.tsv").read_text()
    assert labels == (d / "joined/labels.ko.tsv").read_text() and labels.count("K00102") == 2
    # past the limit a family keeps its commonest labels
    few = pl.DataFrame(
        [(1, "b"), (1, "a"), (1, "b"), (1, "c"), (1, "c"), (1, "c"), (2, "a"), (2, None)],
        orient="row",
        schema=["pfam_accession", "label_ko"],
    )
    bundle.LABELS = 2
    assert bundle.label_names(few, "ko").rows() == [(1, "b", 1), (1, "c", 2), (2, "a", 1)]
    pl.DataFrame(
        [("m", "PF00011", 6, "50M", x[5:55]), ("n", "PF00012", 6, "50M", x[5:55])],
        orient="row",
        schema=["read", "pfam", "hmm_from", "cigar", "residues"],
    ).write_parquet(d / "mx2.parquet")
    got = []
    for b in ("whole", "joined"):
        run(
            "profile",
            "--bundle",
            d / b,
            "--parquet",
            d / "mx2.parquet",
            "--out",
            d / f"q{b}",
            "--hmm",
            d / "mm2.hmm",
        )
        got.append(pl.read_parquet(d / f"q{b}.reads.parquet").sort("read").drop("qid").rows())
    assert got[0] == got[1] and got[0][0][-1] == "MGYP000000000101", got
print("all bundle tests pass")
