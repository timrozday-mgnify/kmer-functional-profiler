"""Extract whole MGnify90 clusters from a MGnify Proteins release into the index's inputs.

One subcommand per pipeline step (``main.nf``):

- ``membership``: clusters with a member in ``--biome`` (substring of the members' biome
  lineages), keeping 1 in ``--sample`` by ``cluster_rep``, -> ``membership.parquet``
  (``cluster_rep``, ``protein_id``), sorted by ``protein_id`` so each range below reads
  only its row groups. Fails if a member lies at or above ``--max-protein-id``.
- ``extract``: sequences (and Pfam hits, unless ``--no-pfam``) of the members in one
  ``protein_id`` range, with their ``cluster_rep`` and ``bucket`` (a hash of
  ``cluster_rep`` modulo ``--buckets``). Fails if any member has no sequence.
- ``merge``: one bucket of the extracted ranges -> ``members.parquet`` (``protein_id``,
  ``cluster_rep``, ``full_length``, ``sequence``, sorted by cluster) and, with ``--pfam``,
  all Pfam hits. Every cluster lies in one bucket, so buckets are independent inputs.
- ``subset``: the clusters of a members table with ``cluster_rep % --sample == 0``; nested,
  since a 1-in-1000 subset is a subset of the 1-in-100 one.
- ``stats``: per-cluster statistics of one bucket for the cost model
  (``kmer_functional_profiler.cost``).
- ``combine``: all buckets' statistics -> ``clusters.parquet``, ``groups.tsv`` (sampled
  k-mers per number of clusters holding them) and ``cost.tsv`` (predicted index size per
  1-in-``--sample`` subset and parameter set).

``--release`` is a local directory or an https prefix (DuckDB reads it remotely).
"""

import argparse
import itertools
import sys
from dataclasses import replace

import duckdb
import polars as pl

from kmer_functional_profiler.cost import STATS_RATE, cluster_stats, predict_cost
from kmer_functional_profiler.index import IndexParams


def connect(args: argparse.Namespace) -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute(f"SET threads = {args.threads}")
    con.execute(f"SET memory_limit = '{args.memory}'")
    con.execute("SET temp_directory = 'duckdb_tmp'")
    con.execute("SET preserve_insertion_order = false")
    return con


def membership(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    con.sql(
        f"""
        SELECT cluster_rep, cluster_member AS protein_id
        FROM read_parquet('{args.release}/mgy_cluster_seqs.parquet')
        WHERE cluster_rep IN (
            SELECT cluster_rep FROM read_parquet('{args.release}/mgy_clusters.parquet')
            WHERE contains(cluster_members_biomes, $biome) AND cluster_rep % $sample = 0
        )
        ORDER BY protein_id
        """,
        params={"biome": args.biome, "sample": args.sample},
    ).write_parquet(args.out)
    top = con.sql(f"SELECT max(protein_id) FROM read_parquet('{args.out}')").fetchone()
    if top and top[0] is not None and top[0] >= args.max_protein_id:
        sys.exit(f"member {top[0]} is at or above --max-protein-id {args.max_protein_id}")


def extract(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    where = "protein_id >= $lo AND protein_id < $hi"
    ids = f"SELECT * FROM read_parquet('{args.membership}') WHERE {where}"
    params = {"lo": args.lo, "hi": args.hi}
    release = args.release
    con.execute(
        f"""
        CREATE TABLE seqs AS
        SELECT m.protein_id, m.cluster_rep, s.full_length, s.sequence,
            hash(m.cluster_rep) % {args.buckets} AS bucket
        FROM ({ids}) m LEFT JOIN (
            SELECT protein_id, full_length, sequence
            FROM read_parquet('{release}/mgy_protein_sequences.parquet') WHERE {where}
        ) s USING (protein_id)
        """,
        params,
    )
    n_members, n_missing = con.sql(
        "SELECT count(*), count(*) - count(sequence) FROM seqs"
    ).fetchone() or (0, 0)
    if n_missing:
        sys.exit(f"{n_missing} of {n_members} members in [{args.lo}, {args.hi}) have no sequence")
    con.sql("SELECT * FROM seqs ORDER BY bucket").write_parquet(f"{args.prefix}.seqs.parquet")
    if args.pfam:
        con.sql(
            "SELECT protein_id, pfam_accession, score, env_from, env_to "
            f"FROM read_parquet('{release}/mgy_proteins_pfam.parquet') "
            f"WHERE {where} AND protein_id IN (SELECT protein_id FROM ({ids}))",
            params=params,
        ).write_parquet(f"{args.prefix}.pfam.parquet")


def merge(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    seqs = f"read_parquet({[f'{p}.seqs.parquet' for p in args.prefixes]})"
    con.sql(
        f"SELECT protein_id, cluster_rep, full_length, sequence FROM {seqs} "
        f"WHERE bucket = {args.bucket} ORDER BY cluster_rep, protein_id"
    ).write_parquet(args.members)
    if args.pfam:
        pfam = f"read_parquet({[f'{p}.pfam.parquet' for p in args.prefixes]})"
        con.sql(f"SELECT * FROM {pfam} ORDER BY protein_id").write_parquet(args.pfam)
    n_members, n_clusters = con.sql(
        f"SELECT count(*), count(DISTINCT cluster_rep) FROM read_parquet('{args.members}')"
    ).fetchone() or (0, 0)
    print(f"bucket {args.bucket}: {n_members} members in {n_clusters} clusters")


def subset(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    members = f"read_parquet('{args.members}')"
    con.sql(
        f"SELECT * FROM {members} WHERE cluster_rep % {args.sample} = 0 "
        "ORDER BY cluster_rep, protein_id"
    ).write_parquet(args.out_members)
    if args.pfam:
        con.sql(
            f"SELECT * FROM read_parquet('{args.pfam}') WHERE protein_id IN "
            f"(SELECT protein_id FROM read_parquet('{args.out_members}')) ORDER BY protein_id"
        ).write_parquet(args.out_pfam)


def stats(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    params = IndexParams(k=args.k, alphabet=args.alphabet)
    clusters, pairs = cluster_stats(args.members, params, args.rate)
    clusters.write_parquet(f"{args.prefix}.clusters.parquet")
    pairs.write_parquet(f"{args.prefix}.pairs.parquet")


def combine(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    pl.scan_parquet([f"{p}.clusters.parquet" for p in args.prefixes]).sink_parquet(
        "clusters.parquet"
    )
    clusters = pl.scan_parquet("clusters.parquet")
    pairs = pl.scan_parquet([f"{p}.pairs.parquet" for p in args.prefixes])
    base = IndexParams(k=args.k, alphabet=args.alphabet)
    groups, rows = [], []
    for sample in args.sample:
        in_sample = pl.col("cluster_rep") % sample == 0
        hist = (
            pairs.filter(in_sample)
            .group_by("hash")
            .len("n_groups")
            .group_by("n_groups")
            .len("hashes")
            .sort("n_groups")
            .collect(engine="streaming")
            .with_columns(sample=pl.lit(sample))
        )
        groups.append(hist)
        pairs_in = hist["n_groups"] * hist["hashes"]
        for t_base, n_min, t_cap, t_dense, max_groups in itertools.product(
            args.t_base, args.n_min, args.t_cap, args.t_dense, args.max_groups
        ):
            params = replace(
                base,
                t_base=t_base,
                n_min=n_min,
                t_cap=t_cap,
                t_dense=t_dense,
                max_groups=max_groups,
            )
            kept = pairs_in.filter(hist["n_groups"] <= max_groups).sum()
            keep = kept / pairs_in.sum() if pairs_in.sum() else 1.0
            cost = predict_cost(clusters.filter(in_sample), params, keep=keep, sets=args.sets)
            rows.append(
                {
                    "sample": sample,
                    "t_base": t_base,
                    "n_min": n_min,
                    "t_cap": t_cap,
                    "t_dense": t_dense,
                    "max_groups": max_groups,
                    "keep": keep,
                    **cost,
                }
            )
    pl.concat(groups).select("sample", "n_groups", "hashes").write_csv("groups.tsv", separator="\t")
    pl.DataFrame(rows).write_csv("cost.tsv", separator="\t", float_precision=6)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--memory", default="4GB")
    sub = parser.add_subparsers(dest="step", required=True)
    p = sub.add_parser("membership")
    p.add_argument("--release", required=True)
    p.add_argument("--biome", default="root")
    p.add_argument("--sample", type=int, default=1)
    p.add_argument("--max-protein-id", type=int, default=11_200_000_000)
    p.add_argument("--out", default="membership.parquet")
    p = sub.add_parser("extract")
    p.add_argument("--release", required=True)
    p.add_argument("--membership", required=True)
    p.add_argument("--lo", type=int, required=True)
    p.add_argument("--hi", type=int, required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--buckets", type=int, default=1)
    p.add_argument("--pfam", action=argparse.BooleanOptionalAction, default=True)
    p = sub.add_parser("merge")
    p.add_argument("--bucket", type=int, default=0)
    p.add_argument("--members", default="members.parquet")
    p.add_argument("--pfam", help="also merge Pfam hits into this file")
    p.add_argument("prefixes", nargs="+")
    p = sub.add_parser("subset")
    p.add_argument("--members", required=True)
    p.add_argument("--pfam")
    p.add_argument("--sample", type=int, required=True)
    p.add_argument("--out-members", required=True)
    p.add_argument("--out-pfam")
    for name in ("stats", "combine"):
        p = sub.add_parser(name)
        p.add_argument("--k", type=int, default=IndexParams.k)
        p.add_argument("--alphabet", default=IndexParams.alphabet)
    stats_p, combine_p = sub.choices["stats"], sub.choices["combine"]
    stats_p.add_argument("--members", required=True)
    stats_p.add_argument("--prefix", required=True)
    stats_p.add_argument("--rate", type=float, default=STATS_RATE)
    combine_p.add_argument("prefixes", nargs="+")
    combine_p.add_argument("--sample", type=int, nargs="+", default=[1])
    for name, kind in (
        ("t_base", float),
        ("n_min", int),
        ("t_cap", float),
        ("t_dense", float),
        ("max_groups", int),
    ):
        default = getattr(IndexParams, name)
        combine_p.add_argument(
            f"--{name.replace('_', '-')}", type=kind, nargs="+", default=[default]
        )
    combine_p.add_argument("--sets", type=float, default=1.0, help="value sets per posting")
    args = parser.parse_args()
    steps = {
        "membership": membership,
        "extract": extract,
        "merge": merge,
        "subset": subset,
        "stats": stats,
        "combine": combine,
    }
    steps[args.step](connect(args), args)


if __name__ == "__main__":
    main()
