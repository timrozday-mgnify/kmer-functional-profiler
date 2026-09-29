"""Extract whole MGnify90 clusters from a MGnify Proteins release into the index's inputs.

One subcommand per pipeline step (``main.nf``):

- ``membership``: clusters with a member in ``--biome`` (substring of the members' biome
  lineages), keeping 1 in ``--sample`` by ``cluster_rep``, -> ``membership/``
  (``cluster_rep``, ``protein_id``), split into ``--shard-width`` ranges of ``protein_id``
  so each range below reads only its own files; a partitioned write, not a sort, since the
  whole release has 5.7 B members. Fails if a member lies at or above ``--max-protein-id``.
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
- ``candidates``, ``bloom``, ``presence``, ``groups``, ``postings``, ``units``,
  ``pack-range``, ``dedup``, ``concat``: the
  partitioned index build over the buckets (``kmer_functional_profiler.partition``); build
  options as ``kmer-functional-profiler index``.
- ``ladder``: nested read subsets of a paired run (``--pairs`` sizes): the pairs ranked
  below n in one shuffle seeded by ``--seed``, so each subset holds the smaller ones.
- ``query-cost``: ``kmer-functional-profiler query --stats`` files named
  ``{index}.{pairs}.{draws}.json`` -> ``query_cost.tsv``, one row per query.

``--release`` is a local directory or an https prefix (DuckDB reads it remotely).
"""

import argparse
import itertools
import json
import subprocess
import sys
from dataclasses import fields, replace
from pathlib import Path

import duckdb
import numpy as np
import polars as pl

from kmer_functional_profiler import partition
from kmer_functional_profiler.cost import STATS_RATE, cluster_stats, predict_cost
from kmer_functional_profiler.index import IndexParams

MEMORY_FRACTION = 0.75  # of the task's memory given to DuckDB


def connect(args: argparse.Namespace) -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute(f"SET threads = {args.threads}")
    # DuckDB's limit leaves out some allocations, and Python needs room too.
    gb = float(args.memory.removesuffix("GB"))
    con.execute(f"SET memory_limit = '{gb * MEMORY_FRACTION:.1f}GB'")
    con.execute("SET temp_directory = 'duckdb_tmp'")
    con.execute("SET preserve_insertion_order = false")
    return con


def membership(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    where = "cluster_rep % $sample = 0"
    if args.biome != "root":  # every cluster is under root: skip the 1.66 B-row join
        where += f"""
            AND cluster_rep IN (
                SELECT cluster_rep FROM read_parquet('{args.release}/mgy_clusters.parquet')
                WHERE contains(cluster_members_biomes, $biome) AND cluster_rep % $sample = 0
            )"""
    params: dict[str, object] = {"sample": args.sample}
    if args.biome != "root":
        params["biome"] = args.biome
    con.execute(
        f"""
        COPY (
            SELECT cluster_rep, cluster_member AS protein_id,
                cluster_member // {args.shard_width} AS shard
            FROM read_parquet('{args.release}/mgy_cluster_seqs.parquet')
            WHERE {where}
        ) TO '{args.out}' (FORMAT parquet, PARTITION_BY (shard), OVERWRITE_OR_IGNORE)
        """,
        params,
    )
    top = con.sql(f"SELECT max(protein_id) FROM {read_membership(args.out)}").fetchone()
    if top and top[0] is not None and top[0] >= args.max_protein_id:
        sys.exit(f"member {top[0]} is at or above --max-protein-id {args.max_protein_id}")


def read_membership(path: str) -> str:
    """Membership files; each holds one ``protein_id`` range, so range filters skip the rest."""
    return f"read_parquet('{path}/**/*.parquet', hive_partitioning = false)"


def extract(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    where = "protein_id >= $lo AND protein_id < $hi"
    ids = f"SELECT cluster_rep, protein_id FROM {read_membership(args.membership)} WHERE {where}"
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


def index_params(args: argparse.Namespace) -> IndexParams:
    return IndexParams(**{f.name: getattr(args, f.name) for f in fields(IndexParams)})


def build(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    params = index_params(args)
    match args.step:
        case "candidates":
            partition.candidates(args.members, args.prefix, params)
        case "bloom":
            partition.bloom(args.prefixes, args.out, args.bits_per_key)
        case "presence":
            partition.presence(
                args.members, args.prefix, args.bloom, args.bucket, args.ranges, params
            )
        case "groups":
            partition.groups(args.paths, args.prefix)
        case "postings":
            partition.postings(args.members, args.prefix, args.groups, params, args.pfam)
        case "units":
            partition.units(args.prefixes, args.bloom, args.ranges, params, args.out, args.pfam)
        case "pack-range":
            partition.pack_range(args.prefixes, args.units, args.range, args.out)
        case "dedup":
            partition.dedup(args.parts, args.range, args.ranges, args.out)
        case "concat":
            partition.concat(args.parts, args.sets, args.units, args.out, params)


def ladder(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    sizes = sorted(args.pairs)
    count = subprocess.Popen(["gzip", "-dc", args.r1], stdout=subprocess.PIPE)
    n_pairs = int(subprocess.run(["wc", "-l"], stdin=count.stdout, capture_output=True).stdout) // 4
    if sizes[-1] > n_pairs:
        sys.exit(f"--pairs {sizes[-1]}: the run has {n_pairs} pairs")
    rank = np.random.default_rng(args.seed).permutation(n_pairs)
    # First subset each pair is in (it is in every larger one too); len(sizes) = none.
    first = np.searchsorted(sizes, rank, side="right").astype(np.uint8).tobytes()
    del rank
    reads = [
        subprocess.Popen(["gzip", "-dc", r], stdout=subprocess.PIPE) for r in (args.r1, args.r2)
    ]
    outs = [
        [
            subprocess.Popen(
                f"gzip -1 > reads.{n}_{mate}.fastq.gz", shell=True, stdin=subprocess.PIPE
            )
            for n in sizes
        ]
        for mate in (1, 2)
    ]
    mates = [zip(*[r.stdout] * 4, strict=False) for r in reads]  # type: ignore[list-item]
    # strict: the mates must have as many reads.
    for i, pair in enumerate(zip(*mates, strict=True)):
        if first[i] < len(sizes):
            r1, r2 = (b"".join(lines) for lines in pair)
            for j in range(first[i], len(sizes)):
                outs[0][j].stdin.write(r1)  # type: ignore[union-attr]
                outs[1][j].stdin.write(r2)  # type: ignore[union-attr]
    for proc in [count, *reads, *outs[0], *outs[1]]:
        if proc.stdin:
            proc.stdin.close()
        if proc.wait():
            sys.exit(f"{proc.args} failed")


def query_cost(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    rows = []
    for path in map(Path, args.paths):
        index, pairs, draws = path.name.removesuffix(".json").split(".")
        stats = json.loads(path.read_text())
        row = {"index": index, "pairs": int(pairs), "draws": int(draws), **stats["counts"]}
        for stage, cost in stats["stages"].items():
            row |= {f"{stage}_{k}": v for k, v in cost.items()}
        rows.append(row)
    table = pl.DataFrame(rows, infer_schema_length=None).sort("index", "pairs", "draws")
    table.write_csv("query_cost.tsv", separator="\t", float_precision=3)


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
    p.add_argument("--shard-width", type=int, default=2**62, help="protein_ids per file")
    p.add_argument("--out", default="membership")
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
    build_steps = (
        "candidates",
        "bloom",
        "presence",
        "groups",
        "postings",
        "units",
        "pack-range",
        "dedup",
        "concat",
    )
    for name in build_steps:
        p = sub.add_parser(name)
        for f in fields(IndexParams):  # the options of `kmer-functional-profiler index`
            flag = f"--{f.name.replace('_', '-')}"
            if isinstance(f.default, bool):
                p.add_argument(flag, action=argparse.BooleanOptionalAction, default=f.default)
            else:
                p.add_argument(flag, type=type(f.default), default=f.default)
    for name in ("candidates", "presence", "postings"):
        sub.choices[name].add_argument("--members", required=True)
        sub.choices[name].add_argument("--prefix", required=True)
    sub.choices["bloom"].add_argument("prefixes", nargs="+")
    sub.choices["bloom"].add_argument("--out", default="bloom")
    sub.choices["bloom"].add_argument("--bits-per-key", type=float, default=10.0)
    sub.choices["presence"].add_argument("--bloom", default="bloom")
    sub.choices["presence"].add_argument("--bucket", type=int, required=True)
    sub.choices["presence"].add_argument("--ranges", type=int, required=True)
    sub.choices["groups"].add_argument("--prefix", required=True)
    sub.choices["groups"].add_argument("paths", nargs="+")
    sub.choices["postings"].add_argument("--pfam")
    sub.choices["postings"].add_argument("groups", nargs="*")
    units_p, range_p, concat_p = (sub.choices[n] for n in ("units", "pack-range", "concat"))
    units_p.add_argument("prefixes", nargs="+")
    units_p.add_argument("--bloom", default="bloom")
    units_p.add_argument("--ranges", type=int, required=True)
    units_p.add_argument("--out", default="units")
    units_p.add_argument("--pfam", action="store_true")
    range_p.add_argument("prefixes", nargs="+")
    range_p.add_argument("--units", default="units")
    range_p.add_argument("--range", type=int, required=True)
    range_p.add_argument("--out", required=True)
    dedup_p = sub.choices["dedup"]
    dedup_p.add_argument("parts", nargs="+", help="pack-range prefixes")
    dedup_p.add_argument("--range", type=int, required=True)
    dedup_p.add_argument("--ranges", type=int, required=True)
    dedup_p.add_argument("--out", required=True)
    concat_p.add_argument("parts", nargs="+", help="pack-range prefixes, in range order")
    concat_p.add_argument("--sets", nargs="+", required=True, help="dedup prefixes, in order")
    concat_p.add_argument("--units", default="units")
    concat_p.add_argument("--out", default="index")
    p = sub.add_parser("ladder")
    p.add_argument("--r1", required=True)
    p.add_argument("--r2", required=True)
    p.add_argument("--pairs", type=int, nargs="+", required=True, help="subset sizes")
    p.add_argument("--seed", type=int, default=1)
    p = sub.add_parser("query-cost")
    p.add_argument("paths", nargs="+")
    args = parser.parse_args()
    steps = {
        "membership": membership,
        "extract": extract,
        "merge": merge,
        "subset": subset,
        "stats": stats,
        "combine": combine,
        **dict.fromkeys(build_steps, build),
        "ladder": ladder,
        "query-cost": query_cost,
    }
    steps[args.step](connect(args), args)


if __name__ == "__main__":
    main()
