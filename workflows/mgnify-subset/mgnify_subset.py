"""Extract whole MGnify90 clusters from a MGnify Proteins release into the index's inputs.

One subcommand per pipeline step (``main.nf``):

- ``membership``: clusters with a member in ``--biome`` (substring of the members' biome
  lineages), keeping 1 in ``--sample`` by ``cluster_rep``, -> ``membership.parquet``
  (``cluster_rep``, ``protein_id``).
- ``extract``: sequences and Pfam hits of those members in one ``protein_id`` range;
  the release files are sorted by ``protein_id``, so each range reads only its row groups.
- ``merge``: join membership with the extracted ranges -> ``members.parquet``
  (``protein_id``, ``cluster_rep``, ``full_length``, ``sequence``) and ``pfam.parquet``.
  Fails if any member has no sequence.

``--release`` is a local directory or an https prefix (DuckDB reads it remotely).
"""

import argparse
import sys

import duckdb


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
        """,
        params={"biome": args.biome, "sample": args.sample},
    ).write_parquet(args.out)


def extract(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    where = "protein_id >= $lo AND protein_id < $hi"
    ids = f"SELECT protein_id FROM read_parquet('{args.membership}') WHERE {where}"
    params = {"lo": args.lo, "hi": args.hi}
    for table, columns, out in (
        ("mgy_protein_sequences", "protein_id, full_length, sequence", "seqs"),
        ("mgy_proteins_pfam", "protein_id, pfam_accession, score, env_from, env_to", "pfam"),
    ):
        con.sql(
            f"SELECT {columns} FROM read_parquet('{args.release}/{table}.parquet') "
            f"WHERE {where} AND protein_id IN ({ids})",
            params=params,
        ).write_parquet(f"{args.prefix}.{out}.parquet")


def merge(con: duckdb.DuckDBPyConnection, args: argparse.Namespace) -> None:
    seqs = f"read_parquet({[f'{p}.seqs.parquet' for p in args.prefixes]})"
    pfam = f"read_parquet({[f'{p}.pfam.parquet' for p in args.prefixes]})"
    con.execute(
        f"""
        CREATE TABLE members AS
        SELECT m.protein_id, m.cluster_rep, s.full_length, s.sequence
        FROM read_parquet('{args.membership}') m LEFT JOIN {seqs} s USING (protein_id)
        """
    )
    n_members, n_missing, n_clusters = con.sql(
        "SELECT count(*), count(*) - count(sequence), count(DISTINCT cluster_rep) FROM members"
    ).fetchone() or (0, 0, 0)
    if n_missing:
        sys.exit(f"{n_missing} of {n_members} members have no sequence; raise --max-protein-id?")
    con.sql("SELECT * FROM members ORDER BY cluster_rep, protein_id").write_parquet(args.members)
    con.sql(f"SELECT * FROM {pfam} ORDER BY protein_id").write_parquet(args.pfam)
    print(f"{n_members} members in {n_clusters} clusters")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--memory", default="4GB")
    sub = parser.add_subparsers(dest="step", required=True)
    p = sub.add_parser("membership")
    p.add_argument("--release", required=True)
    p.add_argument("--biome", default="root")
    p.add_argument("--sample", type=int, default=1)
    p.add_argument("--out", default="membership.parquet")
    p = sub.add_parser("extract")
    p.add_argument("--release", required=True)
    p.add_argument("--membership", required=True)
    p.add_argument("--lo", type=int, required=True)
    p.add_argument("--hi", type=int, required=True)
    p.add_argument("--prefix", required=True)
    p = sub.add_parser("merge")
    p.add_argument("--membership", required=True)
    p.add_argument("--members", default="members.parquet")
    p.add_argument("--pfam", default="pfam.parquet")
    p.add_argument("prefixes", nargs="+")
    args = parser.parse_args()
    {"membership": membership, "extract": extract, "merge": merge}[args.step](connect(args), args)


if __name__ == "__main__":
    main()
