"""Steps of the aai-alpha workflow (``main.nf``): the shape ``alpha`` of the k-mer survival
model that ``aai`` inverts (plan, phase 7, step 23), fitted on MGnify protein pairs.

A pair is a protein P and an MGnify90 cluster C it aligns to. Its survival is P's k-mers
found in the union of C's members' k-mers over ``pin_sum``, C's average member's k-mers:
what the query measures when P's gene is in a sample. Its identity is DIAMOND's, P against
C's representative. ``alpha`` is the least-squares fit of survival against
``query.survival``(identity; alpha) over the pairs.

- ``queries``: ``--n`` random members (seeded) as FASTA in ``--chunks`` files, and every
  cluster's representative as ``reps.faa``.
- ``pairs``: DIAMOND hits (outfmt 6 ``qseqid sseqid pident length qlen slen bitscore``) ->
  one row per (P, C), its best HSP, with both coverages >= ``--min-cov``, P not C's
  representative, at most ``--per-bin`` pairs per identity bin of ``--bin`` from
  ``--min-id``: identities spread evenly, so no band dominates the fit.
- ``survival``: per pair, P's k-mers in C's union and ``pin_sum``, every k-mer (no hash
  threshold: sampling under it is unbiased, all k-mers are more precise). When P is one of
  C's members it is left out of both (a gene new to the index). Clusters are processed in
  batches of about ``--batch-residues`` member residues.
- ``fit``: ``alpha.json`` (``survival``, ``alpha``, ``k``, ``alphabet``, ``fit``), what
  ``kmer-functional-profiler aai-model`` attaches to an index, and ``alpha_strata.tsv``,
  alpha refitted per identity band and per cluster size, to check that one parameter holds.
"""

import argparse
import json
import random
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from scipy.optimize import minimize_scalar

from kmer_functional_profiler import _core
from kmer_functional_profiler.query import survival

HIT_SCHEMA = {
    "protein_id": pl.Int64,
    "cluster_rep": pl.Int64,
    "pident": pl.Float64,
    "length": pl.Int64,
    "qlen": pl.Int64,
    "slen": pl.Int64,
    "bitscore": pl.Float64,
}
MEMBER_COLUMNS = ["protein_id", "cluster_rep", "full_length", "sequence"]
IDENTITY_BANDS = ((0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 1.01))
SIZE_BANDS = ((2, 4), (4, 11), (11, 101), (101, 2**62))  # members, P left out


def _members(paths: list[str]) -> pl.LazyFrame:
    paths = [str(Path(p) / "*.parquet") if Path(p).is_dir() else p for p in paths]
    return pl.scan_parquet(paths).select(MEMBER_COLUMNS)


def _fasta(df: pl.DataFrame, name: str, path: str | Path) -> None:
    df.select(pl.format(">{}\n{}", name, "sequence")).write_csv(
        path, include_header=False, quote_style="never"
    )


def queries(args: argparse.Namespace) -> None:
    members = _members(args.members)
    reps = members.filter(pl.col("protein_id") == pl.col("cluster_rep"))
    reps.select(pl.format(">{}\n{}", "cluster_rep", "sequence")).sink_csv(
        "reps.faa", include_header=False, quote_style="never"
    )
    ids = members.select("protein_id").collect()["protein_id"].to_list()
    chosen = sorted(random.Random(args.seed).sample(ids, min(args.n, len(ids))))
    picked = members.filter(pl.col("protein_id").is_in(chosen)).collect()
    for i, part in enumerate(np.array_split(np.arange(picked.height), args.chunks)):
        _fasta(picked[part], "protein_id", f"queries_{i:03d}.faa")


def pairs(args: argparse.Namespace) -> None:
    hits = pl.concat(
        pl.read_csv(p, separator="\t", has_header=False, schema=HIT_SCHEMA) for p in args.hits
    )
    best = (
        hits.group_by("protein_id", "cluster_rep")
        .agg(pl.all().sort_by("bitscore").last())
        .filter(
            pl.col("protein_id") != pl.col("cluster_rep"),
            pl.col("length") / pl.col("qlen") >= args.min_cov,
            pl.col("length") / pl.col("slen") >= args.min_cov,
            pl.col("pident") / 100 >= args.min_id,
        )
        .select("protein_id", "cluster_rep", identity=pl.col("pident") / 100)
        .with_columns(bin=((pl.col("identity") - args.min_id) / args.bin).floor())
    )
    (
        best.sort("protein_id", "cluster_rep")
        .with_columns(order=pl.struct("protein_id", "cluster_rep").hash(args.seed))  # shuffled
        .filter(pl.col("order").rank("ordinal").over("bin") <= args.per_bin)
        .drop("bin", "order")
        .write_parquet(args.out)
    )


def pair_survival(
    pairs: pl.DataFrame, members: pl.DataFrame, k: int, alphabet: str
) -> pl.DataFrame:
    """``pairs`` (``protein_id`` P, ``cluster_rep`` C, ``identity``) with ``shared`` (P's
    distinct k-mers in C's union), ``pin_sum`` and ``n_members``, P left out of C when it is
    a member; ``members`` holds C's members and every P. ``pin_sum`` averages over the members
    ``p_in`` counts (full-length ones, or all when none is), as the index does.
    ponytail: no P_IN floor for k-mers only fragments hold; it moves pin_sum by < 1 k-mer."""

    def kmers(df: pl.DataFrame) -> pl.DataFrame:
        got = _core.hash_proteins([s.encode() for s in df["sequence"]], k, alphabet=alphabet)
        rows = df.select("protein_id", row=pl.int_range(pl.len(), dtype=pl.UInt64))
        return (
            pl.DataFrame({"row": got["seq"], "hash": got["hash"]})
            .unique()
            .join(rows, on="row")
            .drop("row")
        )

    in_c = members.join(pairs.select("cluster_rep").unique(), on="cluster_rep").with_columns(
        counted=pl.col("full_length") | ~pl.col("full_length").any().over("cluster_rep")
    )
    member_kmers = kmers(in_c).join(in_c.select("protein_id", "cluster_rep"), on="protein_id")
    holders = member_kmers.group_by("cluster_rep", "hash").agg(n=pl.len())
    sizes = (
        member_kmers.group_by("protein_id")
        .agg(n_kmers=pl.len())
        .join(in_c.select("protein_id", "cluster_rep", "counted"), on="protein_id", how="right")
        .with_columns(pl.col("n_kmers").fill_null(0))
    )
    totals = sizes.group_by("cluster_rep").agg(
        total=(pl.col("n_kmers") * pl.col("counted").cast(pl.UInt32)).sum(),
        counted=pl.col("counted").sum(),
        n_members=pl.len(),
    )
    # P's own row in C, if it is a member: left out of the union, the average and the count
    own = pairs.join(sizes, on=["protein_id", "cluster_rep"], how="left").select(
        "protein_id",
        "cluster_rep",
        member=pl.col("counted").is_not_null(),
        own_kmers=(pl.col("n_kmers") * pl.col("counted").cast(pl.UInt32)).fill_null(0),
        own_counted=pl.col("counted").fill_null(False),
    )
    query = members.join(pairs.select("protein_id").unique(), on="protein_id")
    shared = (
        pairs.select("protein_id", "cluster_rep")
        .join(kmers(query), on="protein_id")
        .join(holders, on=["cluster_rep", "hash"])
        .join(own, on=["protein_id", "cluster_rep"])
        .filter(pl.col("n") > pl.col("member").cast(pl.UInt32))
        .group_by("protein_id", "cluster_rep")
        .agg(shared=pl.len())
    )
    return (
        pairs.join(own, on=["protein_id", "cluster_rep"])
        .join(totals, on="cluster_rep")
        .join(shared, on=["protein_id", "cluster_rep"], how="left")
        .select(
            "protein_id",
            "cluster_rep",
            "identity",
            shared=pl.col("shared").fill_null(0),
            pin_sum=(pl.col("total") - pl.col("own_kmers"))
            / (pl.col("counted") - pl.col("own_counted").cast(pl.UInt32)),
            n_members=pl.col("n_members") - pl.col("member").cast(pl.UInt32),
        )
        .filter(pl.col("n_members") > 0, pl.col("pin_sum") > 0)
        .with_columns(survival=pl.col("shared") / pl.col("pin_sum"))
    )


def survival_step(args: argparse.Namespace) -> None:
    pairs_df = pl.read_parquet(args.pairs)
    clusters = pairs_df["cluster_rep"].unique()
    wanted = pl.concat([clusters, pairs_df["protein_id"].unique()]).unique()
    lazy = _members(args.members)
    # every member of C, and P wherever it sits
    members = lazy.filter(
        pl.col("cluster_rep").is_in(clusters.implode())
        | pl.col("protein_id").is_in(wanted.implode())
    ).collect()
    residues = (
        members.filter(pl.col("cluster_rep").is_in(clusters.implode()))
        .group_by("cluster_rep")
        .agg(r=pl.col("sequence").str.len_bytes().sum())
        .sort("cluster_rep")
        .with_columns(batch=pl.col("r").cum_sum() // args.batch_residues)
    )
    parts = []
    for (batch,), group in residues.group_by("batch", maintain_order=True):
        part = pairs_df.join(group.select("cluster_rep"), on="cluster_rep")
        involved = members.filter(
            pl.col("cluster_rep").is_in(group["cluster_rep"].implode())
            | pl.col("protein_id").is_in(part["protein_id"].implode())
        )
        parts.append(pair_survival(part, involved, args.k, args.alphabet))
        print(f"batch {batch}: {part.height} pairs", flush=True)
    pl.concat(parts).write_parquet(args.out)


def fit_alpha(identity: np.ndarray, observed: np.ndarray, k: int) -> tuple[float, float]:
    """(alpha, RMSE) of the least-squares fit of ``observed`` survival against
    survival(``identity``; alpha), alpha in [0.05, 1000] (a log-scale search)."""

    def sse(log_alpha: float) -> float:
        return float(((observed - survival(identity, k, float(np.exp(log_alpha)))) ** 2).sum())

    best: Any = minimize_scalar(sse, bounds=(np.log(0.05), np.log(1000.0)), method="bounded")
    return float(np.exp(best.x)), float(np.sqrt(best.fun / max(len(observed), 1)))


def fit(args: argparse.Namespace) -> None:
    pairs_df = pl.concat(pl.read_parquet(p) for p in args.survival)
    strata = [("all", pl.lit(True))]
    strata += [
        (f"identity {lo}-{min(hi, 1.0)}", pl.col("identity").is_between(lo, hi, closed="left"))
        for lo, hi in IDENTITY_BANDS
    ]
    strata += [
        (
            f"members {lo}-{hi - 1}" if hi < 2**62 else f"members {lo}+",
            pl.col("n_members").is_between(lo, hi, closed="left"),
        )
        for lo, hi in SIZE_BANDS
    ]
    rows = []
    for name, where in strata:
        part = pairs_df.filter(where)
        if name == "all" or part.height >= args.min_stratum:
            a, s = part["identity"].to_numpy(), part["survival"].to_numpy()
            alpha, rmse = fit_alpha(a, s, args.k)
            independent = float(np.sqrt(((s - a**args.k) ** 2).mean()))
            rows.append({"stratum": name, "n": part.height, "alpha": alpha, "rmse": rmse,
                         "rmse_independent": independent})  # fmt: skip
    pl.DataFrame(rows).write_csv(args.strata_out, separator="\t")
    model = {
        "survival": "regional_gamma",
        "alpha": rows[0]["alpha"],
        "k": args.k,
        "alphabet": args.alphabet,
        "fit": {key: rows[0][key] for key in ("n", "rmse", "rmse_independent")},
    }
    Path(args.out).write_text(json.dumps(model, indent=1) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    sub = parser.add_subparsers(dest="step", required=True)
    p = sub.add_parser("queries")
    p.add_argument("--members", required=True, nargs="+", help="parquet files or directories")
    p.add_argument("--n", type=int, default=200_000)
    p.add_argument("--chunks", type=int, default=20)
    p.add_argument("--seed", type=int, default=1)
    p = sub.add_parser("pairs")
    p.add_argument("--hits", required=True, nargs="+")
    p.add_argument("--min-cov", type=float, default=0.8)
    p.add_argument("--min-id", type=float, default=0.6)
    p.add_argument("--bin", type=float, default=0.02)
    p.add_argument("--per-bin", type=int, default=5000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", default="pairs.parquet")
    p = sub.add_parser("survival")
    p.add_argument("--pairs", required=True)
    p.add_argument("--members", required=True, nargs="+")
    p.add_argument("--k", type=int, default=11)
    p.add_argument("--alphabet", default="protein")
    p.add_argument("--batch-residues", type=int, default=50_000_000)
    p.add_argument("--out", default="survival.parquet")
    p = sub.add_parser("fit")
    p.add_argument("--survival", required=True, nargs="+")
    p.add_argument("--k", type=int, default=11)
    p.add_argument("--alphabet", default="protein")
    p.add_argument("--min-stratum", type=int, default=200)
    p.add_argument("--out", default="alpha.json")
    p.add_argument("--strata-out", default="alpha_strata.tsv")
    args = parser.parse_args()
    {"queries": queries, "pairs": pairs, "survival": survival_step, "fit": fit}[args.step](args)


if __name__ == "__main__":
    main()
