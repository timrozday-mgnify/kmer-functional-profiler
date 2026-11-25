"""Is the called cluster a best one when the banded score sees only the domain's residues (+ FLANK)?
By cluster, not by
score, as a shorter target scores less. Run: flank.py LABEL (after `FLANK=... lookup scan`)."""

import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze import OUT, queries  # noqa: E402
from zymo import scores  # noqa: E402

sc = scores("bin")
c = pl.read_csv(OUT / "scan.tsv", separator="\t", schema_overrides={"units": pl.String}).join(
    queries().select("qid", "set", "run", "overall"), on="qid"
)
rows = []
for qid, s, run, best, units in c.select("qid", "set", "run", "overall", "units").iter_rows():
    u, x = sc(qid)
    called = (
        np.array([int(i) for i in units.split(",")], np.uint32) if units else np.zeros(0, np.uint32)
    )
    k = np.searchsorted(u, called)
    v = x[k] if len(k) else np.zeros(1, np.int16)
    rows.append(
        (
            s,
            run in ("ERR4334731", "ERR5866877"),
            bool(v.max() == best),
            bool(v.max() >= 0.95 * best),
            float((v == best).mean()),
            len(k),
        )
    )
d = pl.DataFrame(
    rows, schema=["set", "hard", "any_best", "near", "share_best", "ties"], orient="row"
)
for s, n, a, m, b, t in (
    d.group_by("set")
    .agg(
        pl.len(),
        pl.col("any_best").mean(),
        pl.col("near").mean(),
        pl.col("share_best").mean(),
        pl.col("ties").median(),
    )
    .sort("set")
    .iter_rows()
):
    print(
        sys.argv[1],
        s,
        n,
        f"a best cluster called {a:.4f}, within 5% {m:.4f}, share of the called that are best "
        f"{b:.3f}, median called {t:.0f}",
    )
