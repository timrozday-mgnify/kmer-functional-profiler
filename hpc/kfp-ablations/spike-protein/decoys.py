"""Reversed decoys on spike-protein (plan, phase 13): per sample, alignments to decoys, mates
whose best alignment is a decoy (and of those, mates that also align to a target: the reads
decoys take from targets), and the EM's pairs by source kind. Run from this directory.
"""

import glob

import polars as pl

for f in sorted(glob.glob("results-*/*/pileup/assignments.parquet")):
    d = f.rsplit("/", 1)[0]
    a = pl.scan_parquet(f).select("read_name", "mate", "role", "score")
    best = a.sort("score", descending=True).group_by("read_name", "mate").first()
    n = (
        a.select(alignments=pl.len(), decoy=(pl.col("role") == "decoy").sum())
        .collect()
        .row(0, named=True)
    )
    db = best.filter(pl.col("role") == "decoy").select("read_name", "mate")
    tgt = a.filter(pl.col("role") == "target").select("read_name", "mate").unique()
    n["decoy_best"] = db.select(pl.len()).collect().item()
    n["decoy_best_also_on_target"] = (
        db.join(tgt, on=["read_name", "mate"]).select(pl.len()).collect().item()
    )
    em = (
        pl.read_csv(d + "/sources.tsv", separator="\t")
        .group_by("kind")
        .agg(pl.col("pairs").sum().round(1))
    )
    print(d, n, dict(em.iter_rows()))
