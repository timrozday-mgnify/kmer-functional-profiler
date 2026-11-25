"""Tables of the anchored-lookup comparison, from the Rust prototype's output -> OUT/report/*.tsv.

Similarity of a read to its best target is score / self score (BLOSUM62, Smith-Waterman): about its
identity, and 1 for an exact copy.
"""

import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from prep import AA, OUT  # noqa: E402

REP = OUT / "report"
BINS = [0.0, 0.4, 0.6, 0.8, 0.95]
LABELS = ["<.4", ".4-.6", ".6-.8", ".8-.95", "≥.95"]
pl.Config.set_tbl_rows(300)
pl.Config.set_tbl_cols(30)
pl.Config.set_tbl_width_chars(220)
pl.Config.set_fmt_str_lengths(30)


def blosum() -> np.ndarray:
    from Bio.Align import substitution_matrices

    m = substitution_matrices.load("BLOSUM62")
    return np.array([[m[x][y] for y in AA] for x in AA])


def queries() -> pl.DataFrame:
    q = pl.read_parquet(OUT / "queries.parquet").sort("qid")
    d = np.diag(blosum())
    code = {a: i for i, a in enumerate(AA)}
    self_scores = [int(sum(d[code[c]] for c in s if c in code)) for s in q["sequence"]]
    best = np.fromfile(OUT / "bin/truth.best", dtype=np.dtype([("s", "<i2"), ("n", "<u4")]))
    q = q.with_columns(
        self_score=pl.Series(self_scores),
        overall=pl.Series(best["s"]),
        overall_ties=pl.Series(best["n"]),
    ).with_columns(frac=pl.col("overall") / pl.col("self_score"))
    return q.with_columns(
        bin=pl.col("frac").cut(BINS[1:], labels=LABELS, left_closed=True).cast(str),
        qid=pl.col("qid").cast(pl.Int64),
    ).select(
        "qid",
        "set",
        "run",
        "identity",
        "pfam_accession",
        "source_unit",
        "self_score",
        "overall",
        "overall_ties",
        "frac",
        "bin",
        length=pl.col("sequence").str.len_chars(),
    )


def cand(run: str, q: pl.DataFrame) -> pl.DataFrame:
    c = pl.read_csv(OUT / run / "cand.tsv", separator="\t")
    info = pl.read_csv(OUT / run / "queries.tsv", separator="\t")
    return c.join(q, on="qid").join(info, on="qid")


def summary(c: pl.DataFrame, by: list[str]) -> pl.DataFrame:
    return c.group_by("option", *by).agg(
        n=pl.len(),
        best=pl.col("best").mean(),
        near=pl.col("near").mean(),
        tied=pl.col("tied").mean(),
        cands_median=pl.col("cands").median(),
        cands_mean=pl.col("cands").mean(),
        hits_mean=pl.col("hits").mean(),
    )


def reach(q: pl.DataFrame) -> pl.DataFrame:
    """How often the read's best target over all clusters is one its family's rows can reach."""
    rows = []
    for run, name in (
        ("full_reps_ga", "reps, GA"),
        ("full_reps", "reps, GA + weak"),
        ("full_members", "members, GA"),
    ):
        info = pl.read_csv(OUT / run / "queries.tsv", separator="\t")
        j = q.join(info, on="qid", how="left").with_columns(pl.col("best_score").fill_null(0))
        j = j.with_columns(
            same=pl.col("best_score") == pl.col("overall"),
            near=pl.col("best_score") >= 0.95 * pl.col("overall"),
            none=pl.col("best_score") < 0.8 * pl.col("overall"),
        )
        rows.append(
            j.group_by("set", "bin")
            .agg(pl.len(), pl.col("same").mean(), pl.col("near").mean(), pl.col("none").mean())
            .with_columns(arm=pl.lit(name))
        )
        rows.append(
            j.group_by("set")
            .agg(pl.len(), pl.col("same").mean(), pl.col("near").mean(), pl.col("none").mean())
            .with_columns(arm=pl.lit(name), bin=pl.lit("all"))
            .select(rows[-1].columns)
        )
    return pl.concat(rows).sort("set", "bin", "arm")


def index(run: str) -> pl.DataFrame:
    i = pl.read_csv(OUT / run / "index.tsv", separator="\t")
    return (
        i.group_by("option")
        .agg(pl.col("units", "rows", "residues", "postings", "keys").sum())
        .with_columns(
            postings_per_row=pl.col("postings") / pl.col("rows"),
            postings_per_residue=pl.col("postings") / pl.col("residues"),
            units_per_key=pl.col("postings") / pl.col("keys"),
        )
        .sort("option")
    )


def scaling() -> pl.DataFrame:
    """hits and candidates ~ a N^g per option, over families and unit subsamples (1, 1/4, 1/16,
    1/64)."""
    parts = []
    for sub in (0, 2, 4, 6):
        run = "full_reps" if sub == 0 else f"sub{sub}_reps"
        info = pl.read_csv(OUT / run / "queries.tsv", separator="\t").select(
            "qid", "family_units", "best_score", "self_bits"
        )
        q = queries().select("qid", "pfam_accession", "set")
        c = (
            pl.read_csv(OUT / run / "cand.tsv", separator="\t")
            .join(info, on="qid")
            .join(q, on="qid")
        )
        parts.append(
            c.group_by("option", "pfam_accession")
            .agg(
                n=pl.len(),
                N=pl.col("family_units").first(),
                hits=pl.col("hits").mean(),
                cands=pl.col("cands").mean(),
                best=pl.col("best").mean(),
            )
            .with_columns(sub=pl.lit(sub))
        )
    return pl.concat(parts)


def main() -> None:
    REP.mkdir(exist_ok=True)
    q = queries()
    q.write_parquet(REP / "queries.parquet")
    print(
        q.group_by("set", "bin").agg(pl.len(), pl.col("overall_ties").median()).sort("set", "bin")
    )
    r = reach(q)
    r.write_csv(REP / "reach.tsv", separator="\t")
    print(r)
    for run in ("full_reps", "full_reps_ga", "full_members"):
        if not (OUT / run / "cand.tsv").exists():
            continue
        c = cand(run, q)
        s = pl.concat(
            [
                summary(c, ["set"]).with_columns(bin=pl.lit("all")),
                summary(c, ["set", "bin"]).select(*summary(c, ["set"]).columns, "bin"),
            ]
        )
        s.write_csv(REP / f"cand_{run}.tsv", separator="\t")
        index(run).write_csv(REP / f"index_{run}.tsv", separator="\t")
        if run == "full_reps":
            w = s.filter(pl.col("bin") == "all").pivot(
                on="set",
                index="option",
                values=["best", "tied", "cands_median", "cands_mean", "hits_mean"],
            )
            print(w.sort("option"))
            print(index(run))


if __name__ == "__main__":
    main()
