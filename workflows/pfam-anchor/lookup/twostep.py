"""The members pass end to end (scan MEMBERS=1) against representatives only -> report/twostep.tsv.

Simulated reads (windows of members, their own protein left out): is the source cluster among the
called ties? Real reads have no members truth: does the call reach the exhaustive best over
representatives, and how often does the called cluster change? Run after:
  THREADS=4 FLANK=10 MEMBERS=0 lookup scan bin members2/scan_m0.tsv 2000:1000,500:200 0
  THREADS=4 FLANK=10 MEMBERS=1 lookup scan bin members2/scan_m1.tsv 2000:1000,500:200 0
"""

import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze import OUT, REP, queries  # noqa: E402

D = OUT / "members2"


def calls(m: int) -> pl.DataFrame:
    return pl.read_csv(
        D / f"scan_m{m}.tsv", separator="\t", schema_overrides={"units": pl.String}
    ).with_columns(
        arm=pl.lit("reps + members" if m else "reps"), units=pl.col("units").str.split(",")
    )


def main() -> None:
    src = np.fromfile(OUT / "bin/q.src", "<u4")
    q = (
        queries()
        .select("qid", "set", "identity", "overall")
        .with_columns(src=pl.Series(src).cast(pl.String))
    )
    c = pl.concat([calls(0), calls(1)]).join(q, on="qid")
    c = c.with_columns(
        source=pl.col("units").list.contains(pl.col("src")),
        reaches=pl.col("score") >= pl.col("overall"),
    )
    by = ["arm", "k", "second"]
    sim = c.filter(pl.col("set") == "sim")
    real = c.filter(pl.col("set") == "real")
    t = (
        sim.group_by(by)
        .agg(
            pl.len().alias("sim_reads"),
            pl.col("source").mean().alias("sim_source_in_ties"),
            pl.col("ties").median().alias("sim_ties"),
        )
        .join(
            real.group_by(by).agg(
                pl.len().alias("real_reads"),
                pl.col("reaches").mean().alias("real_reaches_reps_best"),
                pl.col("ties").median().alias("real_ties"),
            ),
            on=by,
        )
        .sort("k", "arm")
    )
    # real reads whose called clusters change with members
    a, b = (
        real.filter(pl.col("arm") == x).select("k", "qid", pl.col("units").list.sort())
        for x in ("reps", "reps + members")
    )
    changed = (
        a.join(b, on=["k", "qid"])
        .group_by("k")
        .agg((pl.col("units") != pl.col("units_right")).mean().alias("real_call_changed"))
    )
    t = t.join(changed, on="k")
    t.write_csv(REP / "twostep.tsv", separator="\t")
    pl.Config.set_tbl_cols(12)
    print(t)
    print(
        sim.filter(pl.col("k") == 2000)
        .pivot(on="arm", index="identity", values="source", aggregate_function="mean")
        .sort("identity", descending=True)
    )


if __name__ == "__main__":
    main()
