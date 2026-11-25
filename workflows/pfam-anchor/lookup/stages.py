"""Follow-up stages: how many candidates each stage must keep to keep the best target ->
report/stages.tsv."""

import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze import OUT, REP, queries  # noqa: E402

q = queries().select("qid", "set", "bin", "overall")
f = pl.read_csv(
    OUT / "follow/follow.tsv", separator="\t", schema_overrides={"option": pl.String}
).join(q, on="qid")
s = (
    f.group_by("option", "stage", "rule", "set")
    .agg(
        n=pl.len(),
        rows=pl.col("rows").mean(),
        kept=pl.col("kept").mean(),
        kept_median=pl.col("kept").median(),
        best=pl.col("best").mean(),
        near=pl.col("near").mean(),
    )
    .sort("option", "stage", "rule", "set")
)
s.write_csv(REP / "stages.tsv", separator="\t")
for o in ("bits8/h1", "top500", "1100101/h1"):
    for r in s.filter((pl.col("option") == o) & (pl.col("set") == "real")).iter_rows(named=True):
        print(
            f"{o:12s} {r['stage']:9s} {r['rule']:8s} best {r['best']:.4f} near {r['near']:.4f} "
            f"kept {r['kept']:.1f} med {r['kept_median']:.0f} of {r['rows']:.0f}"
        )
c = pl.read_csv(
    OUT / "follow/calls.tsv",
    separator="\t",
    schema_overrides={"option": pl.String, "units": pl.String},
).join(q, on="qid")
print(
    c.group_by("option", "set")
    .agg(
        pl.len(),
        right=(pl.col("score") == pl.col("overall")).mean(),
        near=(pl.col("score") >= 0.95 * pl.col("overall")).mean(),
        gapped_equal=(pl.col("gapped_equal") > 0).mean(),
        diag_equal=(pl.col("diag_equal") > 0).mean(),
        ties=pl.col("ties").median(),
    )
    .sort("option", "set")
)
