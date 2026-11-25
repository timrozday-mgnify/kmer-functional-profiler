"""The bundle against the prototype scan: same rows, no collapsing, so the same calls. Fails if any
differ."""

import os
import subprocess
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bundle import LOOKUP  # noqa: E402
from prep import OUT, SLICE  # noqa: E402

B = OUT / "bundle_check"
subprocess.run(["rm", "-rf", B], check=True)
B.mkdir()
m = pl.read_parquet(
    SLICE / "members.parquet", columns=["protein_id", "cluster_rep", "sequence"]
).filter(pl.col("protein_id") == pl.col("cluster_rep"))
m.select("protein_id", "sequence").write_parquet(B / "seqs.parquet")
pl.read_parquet(OUT / "domains.parquet").join(
    m.select("protein_id"), on="protein_id", how="semi"
).write_parquet(B / "hits.parquet")
q = pl.read_parquet(OUT / "queries.parquet").sort("qid")


# the prototype's queries in the pipeline's schema: a simulated read's CIGAR starts inside its
# window, so give it as an insert
def cg(af: int, c: str, n: int) -> str:
    import re

    used = sum(int(k) for k, o in re.findall(r"(\d+)([MID])", c) if o != "D")
    return (
        (f"{af - 1}I" if af > 1 else "")
        + c
        + (f"{n - af + 1 - used}I" if n - af + 1 - used else "")
    )


q = q.with_columns(
    cigar=pl.Series(
        [
            cg(*r)
            for r in q.select("ali_from", "cigar", pl.col("sequence").str.len_chars()).iter_rows()
        ]
    )
)
q.select(
    read=pl.col("qid").cast(pl.String),
    pfam=pl.format("PF{}", pl.col("pfam_accession").cast(pl.String).str.zfill(5)),
    hmm_from="hmm_from",
    cigar="cigar",
    residues="sequence",
).write_parquet(B / "sample.parquet")
py = [sys.executable, Path(__file__).parent / "bundle.py"]
subprocess.run(
    [
        *py,
        "build",
        "--hits",
        B / "hits.parquet",
        "--seqs",
        B / "seqs.parquet",
        "--hmm",
        OUT / "hmms.hmm",
        "--out",
        B / "b",
        "--collapse",
        "0",
    ],
    check=True,
)
subprocess.run([*py, "verify", B / "b"], check=True)
subprocess.run(
    [
        *py,
        "profile",
        "--bundle",
        B / "b",
        "--parquet",
        B / "sample.parquet",
        "--out",
        B / "p",
        "--hmm",
        OUT / "hmms.hmm",
    ],
    check=True,
)
subprocess.run(
    [LOOKUP, "scan", OUT / "bin", B / "scan.tsv", "500:200"],
    check=True,
    env={**os.environ, "FLANK": "10", "THREADS": "4"},
    stderr=subprocess.DEVNULL,
)
s = pl.read_csv(B / "scan.tsv", separator="\t", schema_overrides={"units": pl.String}).select(
    qid=pl.col("qid").cast(pl.String), score0="score", ties0="ties", units="units"
)
units = pl.read_parquet(OUT / "bin/units.parquet")["unit"].to_list()
r = pl.read_parquet(B / "p.reads.parquet").join(s, left_on="read", right_on="qid")
assert r.height == q.height, (r.height, q.height)
bad = r.filter(pl.col("score") != pl.col("score0"))
named = [
    int(p[4:]) in {units[int(u)] for u in us.split(",")} or t0 > 200
    for p, us, t0 in r.select("protein", "units", "ties0").iter_rows()
    if p
]
print(
    f"{r.height} reads: {bad.height} scores differ; {named.count(False)} named proteins are not "
    "among the scan's best"
)
if bad.height or not all(named):
    print(bad.head())
    sys.exit(1)
