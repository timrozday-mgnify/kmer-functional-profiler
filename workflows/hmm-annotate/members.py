"""MGnify90 members' candidate models, from their representatives' hits.

A member is at least 90% identical to its representative, so the models worth HMMER's time on it are
those its representative nearly passes (``annotate.py search --near``).

- ``pairs``: a range of representatives' near hits -> their members' (protein_id, model), in
  ``--out``/shard=S/ by the member's id (S = protein_id // ``--width``), and ``--shards``, the S.
  With ``--check-sample`` 1 in that many members is there without a candidate too (model null);
  with ``--all`` every one is, and every protein of the clusters under ``--min-cluster``, whose
  representatives had no search of their own.
- ``fasta``: one shard of every batch's pairs -> ``--out`` FASTA, ``--out``.tsv (its proteins and
  residues) and ``--pairs-out``, the pairs as ``filter join`` streams them, ``--block`` proteins
  at a time.

``--release`` is a local directory or an https prefix (DuckDB reads it remotely).
"""

import argparse
import struct
import sys
from pathlib import Path

import duckdb
import numpy as np


def connect(a: argparse.Namespace) -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute(f"SET threads = {a.threads}")
    # DuckDB's limit leaves out some allocations, and Python needs room too.
    con.execute(f"SET memory_limit = '{float(a.memory.removesuffix('GB')) * 0.75:.1f}GB'")
    con.execute("SET temp_directory = 'duckdb_tmp'")
    return con


def sql_list(files: list[Path]) -> str:
    return "[" + ", ".join(f"'{f}'" for f in files) + "]"


def pairs(a: argparse.Namespace) -> None:
    con = connect(a)
    con.execute("SET preserve_insertion_order = false")
    wanted = "n.model IS NOT NULL"
    if a.all:
        wanted = "true"
    elif a.check_sample:
        wanted += f" OR c.cluster_member % {a.check_sample} = 0"
    reps = "true"
    if a.min_cluster > 1:
        reps = f"""c.cluster_rep IN (
            SELECT cluster_rep FROM read_parquet('{a.release}/mgy_clusters.parquet')
            WHERE cluster_size >= {a.min_cluster} AND cluster_rep BETWEEN {a.lo} AND {a.hi - 1})"""
    a.out.mkdir()
    con.execute(
        f"""
        COPY (
            SELECT c.cluster_member AS protein_id, n.model, c.cluster_member // {a.width} AS shard
            FROM read_parquet('{a.release}/mgy_cluster_seqs.parquet') c
            LEFT JOIN read_parquet({sql_list(a.near)}) n ON n.protein_id = c.cluster_rep
            WHERE c.cluster_rep BETWEEN {a.lo} AND {a.hi - 1} AND c.cluster_rep % {a.sample} = 0
                AND NOT (c.cluster_member = c.cluster_rep AND {reps}) AND ({wanted})
        ) TO '{a.out}' (FORMAT parquet, COMPRESSION zstd, PARTITION_BY (shard), OVERWRITE_OR_IGNORE)
        """
    )
    shards = sorted(int(p.name.removeprefix("shard=")) for p in a.out.glob("shard=*"))
    a.shards.write_text("".join(f"{s}\n" for s in shards))


def fasta(a: argparse.Namespace) -> None:
    files = sorted(f for d in a.pairs for f in (d / f"shard={a.shard}").glob("*.parquet"))
    con = connect(a)
    con.execute(
        f"""
        CREATE TABLE p AS SELECT protein_id, t.number
        FROM read_parquet({sql_list(files)}, hive_partitioning = false) r
        LEFT JOIN (
            SELECT model, row_number() OVER () - 1 AS number
            FROM read_csv('{a.thresholds}', delim='\t', header=true, all_varchar=true)
        ) t USING (model) ORDER BY protein_id
        """
    )
    wanted, known = con.sql("SELECT count(DISTINCT protein_id), count(number) FROM p").fetchone()
    given = con.sql(f"SELECT count(model) FROM read_parquet({sql_list(files)})").fetchone()[0]
    unknown = given - known
    if unknown:
        sys.exit(f"{unknown} pairs' models are not in {a.thresholds}")
    rows = con.cursor().execute(
        f"""
        SELECT protein_id, sequence FROM read_parquet('{a.release}/mgy_protein_sequences.parquet')
        WHERE protein_id BETWEEN {a.shard * a.width} AND {(a.shard + 1) * a.width - 1}
            AND protein_id IN (SELECT protein_id FROM p) ORDER BY protein_id
        """
    )
    n = residues = 0
    with open(a.out, "w") as out, open(a.pairs_out, "wb") as stream:
        while block := rows.fetchmany(a.block):
            out.writelines(f">{i}\n{s}\n" for i, s in block)
            ids = np.array([i for i, _ in block])
            size = sum(len(s) for _, s in block)
            found = con.execute(
                f"SELECT protein_id, number FROM p WHERE number IS NOT NULL "
                f"AND protein_id BETWEEN {ids[0]} AND {ids[-1]} ORDER BY number, protein_id"
            ).fetchnumpy()
            data = np.stack([np.searchsorted(ids, found["protein_id"]), found["number"]], axis=1)
            data = data.astype("<u4").tobytes()
            stream.write(struct.pack("<3Q", len(block), size, len(data)) + data)
            n, residues = n + len(block), residues + size
    if n != wanted:
        sys.exit(f"shard {a.shard}: {n} sequences of {wanted} members")
    Path(f"{a.out}.tsv").write_text(f"proteins\tresidues\n{n}\t{residues}\n")


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--memory", default="4GB")
    sub = p.add_subparsers(required=True)
    x = sub.add_parser("pairs")
    x.set_defaults(run=pairs)
    x.add_argument("--near", type=Path, nargs="+", required=True)
    x.add_argument("--lo", type=int, required=True, help="the representatives' range")
    x.add_argument("--hi", type=int, required=True)
    x.add_argument("--sample", type=int, default=1, help="1 in this many clusters")
    x.add_argument("--check-sample", type=int, default=0)
    x.add_argument("--all", action="store_true")
    x.add_argument("--min-cluster", type=int, default=1)
    x.add_argument("--shards", type=Path, default=Path("shards.txt"))
    f = sub.add_parser("fasta")
    f.set_defaults(run=fasta)
    f.add_argument("--pairs", type=Path, nargs="+", required=True, help="every `pairs` --out")
    f.add_argument("--shard", type=int, required=True)
    f.add_argument("--thresholds", type=Path, required=True)
    f.add_argument("--pairs-out", type=Path, required=True)
    f.add_argument("--block", type=int, default=80_000)
    for s in (x, f):
        s.add_argument("--release", required=True)
        s.add_argument("--width", type=int, required=True, help="protein_ids to a shard")
        s.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.run(a)


if __name__ == "__main__":
    main()
