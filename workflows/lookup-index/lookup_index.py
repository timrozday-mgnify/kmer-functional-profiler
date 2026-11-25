"""lookup-index's tables: hmm-annotate's Pfam domains -> part rows (``bundle.PART``) by family.

- ``groups``: the domains' passing rows per family -> ``groups.tsv`` (model, rows, group), families
  in about ``--groups`` groups of equal rows, a larger family alone; and ``batches.tsv``, the
  annotation's chunks taken ``--batch-chunks`` at a time with their protein_id range. Fails unless
  every chunk of the annotations is there and ``--hmm`` is the file they were made with.
- ``split``: one batch -> ``--out``/group=G/*.parquet. Joins each protein's cluster
  (``--membership``, mgnify-subset's shards) and its labels: ``--label ko=DIR``, an hmm-annotate
  output whose models are the labels; a protein's several are one name, joined by +.
"""

import argparse
import csv
import gzip
import hashlib
import heapq
import sys
from pathlib import Path

import duckdb


def connect(args: argparse.Namespace) -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute(f"SET threads = {args.threads}")
    # DuckDB's limit leaves out some allocations, and Python needs room too.
    con.execute(f"SET memory_limit = '{float(args.memory.removesuffix('GB')) * 0.75:.1f}GB'")
    con.execute("SET temp_directory = 'duckdb_tmp'")
    con.execute("SET preserve_insertion_order = false")
    return con


def chunks(annotation: Path, kind: str) -> list[tuple[int, int, Path]]:
    """(lo, hi, file) of every chunk of an hmm-annotate output, by protein_id."""
    with open(annotation / "ranges.tsv") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    out = [
        (int(r["lo"]), int(r["hi"]), annotation / kind / f"chunk{r['chunk']}.{kind}.parquet")
        for r in rows
    ]
    missing = [p for _, _, p in out if not p.is_file()]
    if missing or not out:
        sys.exit(f"{len(missing)} of {len(out)} chunks of {annotation} are missing: {missing[:3]}")
    return sorted(out)


def sql_list(files: list[Path]) -> str:
    return "[" + ", ".join(f"'{f}'" for f in files) + "]"


def label_dirs(specs: list[str]) -> list[tuple[str, Path]]:
    out = []
    for spec in specs:
        scheme, _, path = spec.partition("=")
        if not scheme.isidentifier() or not path:
            sys.exit(f"--label {spec}: scheme=directory, the scheme letters, digits and _")
        out.append((scheme, Path(path)))
    return out


def check_hmm(hmm: Path, annotation: Path) -> None:
    opener = gzip.open if hmm.suffix == ".gz" else open
    with opener(hmm, "rb") as f:
        digest = hashlib.sha256()
        while data := f.read(1 << 24):
            digest.update(data)
    if digest.hexdigest() != (annotation / "hmm.sha256").read_text().strip():
        sys.exit(f"{hmm} is not the HMM file {annotation} was made with")


def groups(args: argparse.Namespace) -> None:
    check_hmm(args.hmm, args.pfam)
    found = chunks(args.pfam, "domains")
    for _, d in label_dirs(args.label):
        chunks(d, "labels")
    files = sql_list([p for _, _, p in found])
    counts = (
        connect(args)
        .sql(
            f"SELECT model, count(*) n FROM read_parquet({files}) WHERE passes "
            "GROUP BY model ORDER BY n DESC, model"
        )
        .fetchall()
    )
    odd = [m for m, _ in counts if not (m.startswith("PF") and m[2:].isdigit())]
    if odd or not counts:
        sys.exit(f"{args.pfam}: {len(counts)} models, {len(odd)} not Pfam accessions: {odd[:3]}")
    # largest first into the lightest group, which leaves the largest alone if it fills one
    loads = [(0, g) for g in range(min(args.groups, len(counts)))]
    with open(args.out / "groups.tsv", "w") as out:
        out.write("model\trows\tgroup\n")
        for model, n in counts:
            load, g = heapq.heappop(loads)
            out.write(f"{model}\t{n}\t{g}\n")
            heapq.heappush(loads, (load + n, g))
    with open(args.out / "batches.tsv", "w") as out:
        out.write("batch\tlo\thi\n")
        for b, i in enumerate(range(0, len(found), args.batch_chunks)):
            part = found[i : i + args.batch_chunks]
            out.write(f"{b}\t{part[0][0]}\t{part[-1][1]}\n")
    print(f"{sum(n for _, n in counts)} rows, {len(counts)} families, the largest {counts[0]}")


def split(args: argparse.Namespace) -> None:
    con = connect(args)
    files = [p for lo, hi, p in chunks(args.pfam, "domains") if lo < args.hi and hi > args.lo]
    in_range = f"protein_id BETWEEN {args.lo} AND {args.hi - 1}"
    columns, joins = [], []
    for scheme, d in label_dirs(args.label):
        near = [p for lo, hi, p in chunks(d, "labels") if lo < args.hi and hi > args.lo]
        source = (
            f"SELECT protein_id, string_agg(DISTINCT model, '+' ORDER BY model) AS label "
            f"FROM read_parquet({sql_list(near)}) WHERE {in_range} GROUP BY protein_id"
            if near
            else "SELECT 0::BIGINT AS protein_id, NULL::VARCHAR AS label WHERE false"
        )
        columns.append(f", l_{scheme}.label AS label_{scheme}")
        joins.append(f"LEFT JOIN ({source}) l_{scheme} USING (protein_id)")
    rows = f"read_parquet({sql_list(files)})"
    con.execute(
        f"""
        COPY (
            SELECT protein_id, CAST(d.model[3:] AS INTEGER) AS pfam_accession, d.hmm_from, d.cigar,
                d."left", d."right", d.residues, m.cluster_rep{"".join(columns)}, g."group"
            FROM {rows} d
            JOIN read_csv('{args.groups}', delim='\t', header=true) g ON g.model = d.model
            LEFT JOIN (
                SELECT protein_id, cluster_rep
                FROM read_parquet('{args.membership}/shard_*.parquet') WHERE {in_range}
            ) m USING (protein_id)
            {" ".join(joins)}
            WHERE d.passes
        ) TO '{args.out}' (FORMAT parquet, COMPRESSION zstd, PARTITION_BY ("group"))
        """
    )
    expected = con.sql(f"SELECT count(*) FROM {rows} WHERE passes").fetchone()[0]
    if not expected:  # a batch with no domain at all writes nothing
        args.out.mkdir(exist_ok=True)
        return
    n, unclustered = con.sql(
        "SELECT count(*), count(*) FILTER (cluster_rep IS NULL) "
        f"FROM read_parquet('{args.out}/*/*.parquet')"
    ).fetchone()
    if n != expected or unclustered:
        sys.exit(f"{n} rows of {expected} domains, {unclustered} of proteins in no cluster")
    print(f"{n} rows")


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--memory", default="4GB")
    sub = p.add_subparsers(required=True)
    g = sub.add_parser("groups")
    g.set_defaults(run=groups)
    g.add_argument("--hmm", type=Path, required=True)
    g.add_argument("--groups", type=int, default=400)
    g.add_argument("--batch-chunks", type=int, default=50)
    g.add_argument("--out", type=Path, default=Path("."))
    s = sub.add_parser("split")
    s.set_defaults(run=split)
    s.add_argument("--groups", type=Path, required=True, help="groups.tsv")
    s.add_argument("--membership", type=Path, required=True)
    s.add_argument("--lo", type=int, required=True)
    s.add_argument("--hi", type=int, required=True)
    s.add_argument("--out", type=Path, default=Path("rows"))
    for x in (g, s):
        x.add_argument("--pfam", type=Path, required=True, help="hmm-annotate's output, of Pfam")
        x.add_argument("--label", action="append", default=[], metavar="SCHEME=DIR")
    args = p.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
