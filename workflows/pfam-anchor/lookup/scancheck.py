"""The end-to-end scan's calls against the exhaustive truth -> printed."""

import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze import OUT, queries  # noqa: E402

q = queries().select("qid", "set", "run", "bin", "overall")
for k in (500, 100):
    c = (
        pl.read_csv(OUT / "scan.tsv", separator="\t", schema_overrides={"units": pl.String})
        .filter(pl.col("k") == k)
        .join(q, on="qid")
    )
    assert c.height == q.height and (c["score"] <= c["overall"]).all()
    c = c.with_columns(
        right=pl.col("score") == pl.col("overall"), near=pl.col("score") >= 0.95 * pl.col("overall")
    )
    print(
        k,
        c.group_by("set")
        .agg(
            pl.len(),
            (~pl.col("right")).sum().alias("missed"),
            pl.col("right").mean(),
            pl.col("near").mean(),
        )
        .sort("set"),
    )
    h = c.filter(pl.col("run").is_in(["ERR4334731", "ERR5866877"]))
    print("hard libraries", h.height, h["right"].mean())
