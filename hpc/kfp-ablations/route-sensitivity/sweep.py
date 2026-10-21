"""Routing's k and min_kmers against unrouted DIAMOND (plan, phase 13, routing): sweep.tsv,
hist.tsv and decision.tsv.

Per sample, the baseline is a finished unrouted ``--sensitive -F 15`` run (step 3b's
``sensitive_F15`` arm: its ``ref/`` and ``pileup/assignments.parquet``). Routing is computed
at component granularity (each component its own block, the strictest routing), for every
k, and thresholded at every min_kmers; DIAMOND is not rerun. A pair "reaches" a sequence when
it is routed to that sequence's component. Columns:

- ``routed_share``: pairs routed (the rest are dropped); ``searched_share``: residues
  searched over pairs x reference residues, DIAMOND's work against the unrouted run;
- ``best_recall``: pairs with a baseline alignment that reach their best alignment's sequence;
- ``unit_recall``: pairs the baseline EM gives a unit (>= 0.5 on one mate) that reach it;
- ``novel_routed``: pairs the baseline EM puts in the novel class (>= 0.5) still routed;
- with ``--truth`` (spike): ``strain_recall`` and ``sibling_recall``, ``best_recall`` over the
  strains' and the withheld siblings' pairs. A sibling pair not routed is not a loss: no
  target can take it, so it counts as caught (``sibling_caught``); ``sibling_recall`` says
  how many siblings still reach DIAMOND, not how well they are caught.

``hist.tsv``: pairs by the most distinct k-mers they share with one component, per k (for
a per-sample threshold fitted like the query's presence model, if a fixed one does not hold).
``decision.tsv``: the rule (``unit_recall`` >= ``MIN_UNIT_RECALL`` in every sample; fixed by
the user before the run, 2026-10-21), the cheapest passing setting by mean
``searched_share``. The chosen setting's truth metrics come from step 3b's re-run with
routing on.
"""

import argparse
import json
import tempfile
from pathlib import Path

import numpy as np
import polars as pl

from kfp_pileup.reference import REF_FASTA, REFS
from kfp_pileup.route import _records, blocks, kmer_table, pair_hits
from kmer_functional_profiler.index import read_fasta

MIN_UNIT_RECALL = 0.99  # the rule, fixed before the run (user, 2026-10-21)
MATE_SUFFIX = r"/[12]$"


def names(r1: Path) -> pl.Series:
    """Pair names in input order (FastxHits' read ids), the first word, /1 dropped."""
    return pl.Series("read_name", [n.decode() for n, _ in _records(r1)]).str.replace(
        MATE_SUFFIX, ""
    )


def baseline(assignments: Path, seq_block: pl.DataFrame) -> dict[str, pl.DataFrame]:
    """Per pair, the block of its best alignment, of its EM unit and its novel share."""
    a = (
        pl.read_parquet(assignments)
        .with_columns(pl.col("read_name").str.replace(MATE_SUFFIX, ""))
        .join(seq_block, on="seq_id")
    )
    best = a.sort("score", descending=True).group_by("read_name").first()
    unit = (
        a.group_by("read_name", "mate", "unit", "block")
        .agg(pl.col("em_posterior").sum())
        .filter(pl.col("em_posterior") >= 0.5)
        .group_by("read_name")
        .agg(pl.col("block").first())
    )
    novel = a.group_by("read_name").agg(pl.col("em_novel").max()).filter(pl.col("em_novel") >= 0.5)
    return {
        "best": best.select("read_name", "block"),
        "unit": unit,
        "novel": novel.select("read_name"),
    }


def reach(pairs: pl.DataFrame, routed: pl.DataFrame) -> float | None:
    """Share of ``pairs`` (read_name, block) routed to their block."""
    if pairs.height == 0:
        return None
    return pairs.join(routed, on=["read_name", "block"], how="semi").height / pairs.height


def main() -> None:
    parser = argparse.ArgumentParser(description="Routing k and min_kmers against DIAMOND")
    parser.add_argument("--bench", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--ref", type=Path, required=True, help="the baseline's ref/")
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True, help="for the alphabet")
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path, required=True)
    parser.add_argument("--truth", type=Path, help="spike truth.tsv (strain, kind)")
    parser.add_argument("--k", type=int, nargs="+", default=[7, 8, 9, 10, 11])
    parser.add_argument("--min-kmers", type=int, nargs="+", default=[1, 2, 3, 5, 8])
    parser.add_argument("--out", type=Path, default=Path("."))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    alphabet = json.loads((args.index / "meta.json").read_text())["params"]["alphabet"]
    refs = pl.read_csv(args.ref / REFS, separator="\t", infer_schema_length=None)
    sequences = dict(read_fasta(args.ref / REF_FASTA))
    seqs = [sequences[s].encode() for s in refs["seq_id"]]
    block = blocks(refs, 1)  # every component its own block
    seq_block = pl.DataFrame({"seq_id": refs["seq_id"], "block": block})
    letters = (
        pl.DataFrame({"block": block, "length": refs["length"]})
        .group_by("block")
        .agg(letters=pl.col("length").sum())
    )
    total = int(refs["length"].sum())
    pair_names = names(args.r1)
    base = baseline(args.assignments, seq_block)
    groups = {"best": base["best"], "unit": base["unit"]}
    if args.truth:
        kind = pl.read_csv(args.truth, separator="\t").select("strain", "kind")
        labelled = (
            base["best"]
            .with_columns(strain=pl.col("read_name").str.extract(r"^s(\d+)_", 1).cast(pl.Int64))
            .join(kind, on="strain")
        )
        groups |= {
            "strain": labelled.filter(pl.col("kind") == "strain"),
            "sibling": labelled.filter(pl.col("kind") == "sibling"),
        }

    rows, hist = [], []
    for k in args.k:
        with tempfile.TemporaryDirectory(prefix=".table-", dir=args.out) as scratch:
            table = kmer_table(seqs, block, k, alphabet, Path(scratch))
            hits, n_pairs = pair_hits(table, (args.r1, args.r2), k, alphabet)
            del table
        assert n_pairs == pair_names.len(), "pair numbering differs from the FASTQ"
        hits = hits.with_columns(read_name=pair_names.gather(hits["read"].cast(pl.Int64)))
        most = hits.group_by("read").agg(pl.col("kmers").max())
        counts = np.bincount(most["kmers"].to_numpy())
        counts[0] = n_pairs - most.height
        hist += [
            {"bench": args.bench, "sample": args.sample, "k": k, "kmers": i, "pairs": int(c)}
            for i, c in enumerate(counts)
            if c
        ]
        for m in args.min_kmers:
            routed = hits.filter(pl.col("kmers") >= m)
            searched = routed.join(letters, on="block")["letters"].sum()
            row = {
                "bench": args.bench,
                "sample": args.sample,
                "k": k,
                "min_kmers": m,
                "pairs": n_pairs,
                "routed_share": routed["read"].n_unique() / n_pairs,
                "searched_share": searched / (n_pairs * total),
                "novel_routed": base["novel"].join(routed, on="read_name", how="semi").height
                / max(base["novel"].height, 1),
            }
            row |= {f"{g}_recall": reach(p, routed) for g, p in groups.items()}
            rows.append(row)
    tag = f"{args.bench}.{args.sample}"
    pl.DataFrame(rows).write_csv(args.out / f"sweep.{tag}.tsv", separator="\t")
    pl.DataFrame(hist).write_csv(args.out / f"hist.{tag}.tsv", separator="\t")


def decide(out: Path) -> None:
    """decision.tsv over every sample's sweep: the rule in every sample, cheapest first."""
    sweep = pl.concat(
        [pl.read_csv(p, separator="\t") for p in sorted(out.glob("sweep.*.tsv"))],
        how="diagonal_relaxed",
    )
    sweep.write_csv(out / "sweep.tsv", separator="\t")
    truth = [c for c in ("strain_recall", "sibling_recall") if c in sweep.columns]
    decision = (
        sweep.group_by("k", "min_kmers")
        .agg(
            pl.col("unit_recall", "best_recall", *truth).min(),
            pl.col("routed_share", "searched_share", "novel_routed").mean(),
        )
        .with_columns(passes=pl.col("unit_recall") >= MIN_UNIT_RECALL)
        .sort("searched_share")
    )
    decision.write_csv(out / "decision.tsv", separator="\t")
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=200):
        print(decision)
    chosen = decision.filter("passes")
    print("chosen:", chosen.select("k", "min_kmers").row(0) if chosen.height else "none")


if __name__ == "__main__":
    import sys

    if sys.argv[1:2] == ["decide"]:
        decide(Path(sys.argv[2]) if len(sys.argv) > 2 else Path("."))
    else:
        main()
