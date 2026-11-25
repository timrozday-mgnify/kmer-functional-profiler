"""What collapsing the targets costs a read: similarity to the best kept target ->
report/annot_collapse.tsv."""

import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1]))
from analyze import OUT, REP, queries  # noqa: E402
from zymo import family_units, scores  # noqa: E402

q = queries().select("qid", "set", "run", "self_score", "overall", "pfam_accession")
units, accs = family_units()
sc = scores("bin")
c = (
    pl.read_csv(OUT / "bin/collapse.tsv", separator="\t")
    .group_by("t")
    .agg(pl.col("rows", "kept", "residues", "residues_kept").sum())
)
out = []
hard = ("ERR4334731", "ERR5866877")
for t in (0.9, 0.8, 0.7, 0.6, 0.5):
    raw, at, masks = np.fromfile(OUT / f"bin/collapse.{t:g}.units", "<u4"), 0, []
    for u in units:
        n = int(raw[at])
        masks.append(np.isin(u, raw[at + 1 : at + 1 + n]))
        at += 1 + n
    rows = []
    for qid, s, run, self_score, best, acc in q.iter_rows():
        _, x = sc(qid)
        rows.append(
            (s, run in hard, best / self_score, float(x[masks[accs.index(acc)]].max()) / self_score)
        )
    d = pl.DataFrame(rows, schema=["set", "hard", "full", "collapsed"], orient="row")
    for name, e in (
        ("real", pl.col("set") == "real"),
        ("real, hard libraries", (pl.col("set") == "real") & pl.col("hard")),
    ):
        x = d.filter(e)
        out.append(
            (
                t,
                name,
                x["full"].mean(),
                x["collapsed"].mean(),
                (x["full"] >= 0.6).mean(),
                (x["collapsed"] >= 0.6).mean(),
                (x["collapsed"] >= 0.4).mean(),
            )
        )
o = pl.DataFrame(
    out,
    schema=[
        "t",
        "reads",
        "similarity_full",
        "similarity_collapsed",
        "full_over_.6",
        "collapsed_over_.6",
        "collapsed_over_.4",
    ],
    orient="row",
)
o = o.join(
    c.with_columns(
        rows_left=pl.col("kept") / pl.col("rows"),
        residues_left=pl.col("residues_kept") / pl.col("residues"),
    ).select("t", "rows_left", "residues_left"),
    on="t",
)
o.write_csv(REP / "annot_collapse.tsv", separator="\t")
pl.Config.set_tbl_cols(12)
pl.Config.set_tbl_width_chars(200)
print(o)
