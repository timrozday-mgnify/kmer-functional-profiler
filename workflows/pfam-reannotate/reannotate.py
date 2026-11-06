"""Re-annotate MGnify proteins with Pfam-A, keeping matches below the gathering thresholds.

One subcommand per pipeline step (``main.nf``). The release's own Pfam table holds only hits
at GA, so this is where the weak ones come from.

- ``ga``: an HMM file's accessions and GA thresholds -> ``ga.tsv``.
- ``plan``: the targets (``reps``: 90% cluster representatives; ``all``: every protein; 1 in
  ``--sample`` of them by id) -> ``ranges.tsv``, ``protein_id`` ranges of about
  ``--chunk-proteins`` each.
- ``fasta``: one range's sequences -> ``--out`` FASTA with ``protein_id`` as the name, and
  ``--out``.tsv with its proteins and residues. Fails if the count differs from the plan's.
  ``--reverse`` reverses each sequence: the null arm (same lengths and composition).
- ``parse``: one ``hmmsearch --domtblout`` and its ``--notextw`` text output -> hits Parquet (one
  row per domain; the release table's column names, plus ``domain``, ``seq_score``,
  ``seq_evalue``, ``c_evalue``, ``ali_from``, ``ali_to``, ``cigar``, ``above_ga``, and ``pp`` with
  ``--pp``). ``cigar`` is over the alignment (M match, I residue in the protein only, D model
  position only): with ``ali_from`` and the searched sequence it rebuilds the aligned protein,
  and with ``hmm_from`` and the HMM the aligned model. ``pp`` is HMMER's posterior probability
  line, a character per alignment column. Also writes a
  one-row statistics TSV. ``above_ga`` is what ``--cut_ga`` would report, except that domtblout
  rounds scores to 0.1 bit, so a domain just under GA can read as above it (8 of 3,261 on 5,000
  MGnify proteins; ``scores`` equal to the threshold).
- ``summary``: the statistics of every chunk -> one total row per arm.
- ``compare``: both arms' hits -> ``compare.tsv``, domains per 5-bit score bin and side of GA.

``--release`` is a local directory or an https prefix (DuckDB reads it remotely).
"""

import argparse
import csv
import gzip
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Final

import duckdb
import polars as pl

BIN_WIDTH: Final = 100_000  # protein_ids per planning bin: at most this many targets in one
BATCH: Final = 100_000  # rows fetched at a time when writing FASTA
SCORE_BIN: Final = 5  # bits per bin of ``compare``
SCHEMA: Final = {
    "protein_id": pl.Int64,
    "pfam_accession": pl.Int32,
    "seq_evalue": pl.Float64,
    "seq_score": pl.Float32,
    "c_evalue": pl.Float64,
    "i_evalue": pl.Float64,
    "score": pl.Float32,
    "hmm_from": pl.Int32,
    "hmm_to": pl.Int32,
    "env_from": pl.Int32,
    "env_to": pl.Int32,
    "domain": pl.Int16,
    "ali_from": pl.Int32,
    "ali_to": pl.Int32,
}
DOMTBL_COLUMNS: Final = (0, 4, 6, 7, 11, 12, 13, 15, 16, 19, 20, 9, 17, 18)  # SCHEMA's, in order


def connect(args: argparse.Namespace) -> duckdb.DuckDBPyConnection:
    con = duckdb.connect()
    con.execute(f"SET threads = {args.threads}")
    # DuckDB's limit leaves out some allocations, and Python needs room too.
    con.execute(f"SET memory_limit = '{float(args.memory.removesuffix('GB')) * 0.75:.1f}GB'")
    con.execute("SET temp_directory = 'duckdb_tmp'")
    return con


def accession(acc: str) -> int:
    """``PF00001.31`` -> 1, as the release table stores it."""
    return int(acc.removeprefix("PF").split(".")[0])


def targets(release: str, kind: str) -> tuple[str, str]:
    """The id column's table and name for one kind of target."""
    if kind == "reps":
        return f"read_parquet('{release}/mgy_clusters.parquet')", "cluster_rep"
    return f"read_parquet('{release}/mgy_cluster_seqs.parquet')", "cluster_member"


def ga(args: argparse.Namespace) -> None:
    acc = None
    with open(args.hmm) as f, open(args.out, "w") as out:
        out.write("pfam_accession\tga_seq\tga_dom\n")
        for line in f:
            if line.startswith("ACC"):
                acc = accession(line.split()[1])
            elif line.startswith("GA"):
                seq, dom = line.split()[1:3]
                out.write(f"{acc}\t{seq}\t{dom.rstrip(';')}\n")


def plan(args: argparse.Namespace) -> None:
    table, column = targets(args.release, args.targets)
    bins = (
        connect(args)
        .sql(
            f"SELECT {column} // {args.bin_width} AS b, count(*) FROM {table} "
            f"WHERE {column} % {args.sample} = 0 GROUP BY b ORDER BY b"
        )
        .fetchall()
    )
    with open(args.out, "w") as out:
        out.write("chunk\tlo\thi\tproteins\n")
        chunk, lo, n = 0, None, 0
        for b, count in bins:
            lo = b * args.bin_width if lo is None else lo
            n += count
            if n >= args.chunk_proteins or b == bins[-1][0]:
                out.write(f"{chunk}\t{lo}\t{(b + 1) * args.bin_width}\t{n}\n")
                chunk, lo, n = chunk + 1, None, 0


def fasta(args: argparse.Namespace) -> None:
    # protein_id sorted in the release's Parquet: a range reads only its own row groups.
    in_range = f"BETWEEN {args.lo} AND {args.hi - 1}"
    where = f"protein_id {in_range} AND protein_id % {args.sample} = 0"
    if args.targets == "reps":
        table, column = targets(args.release, "reps")
        where += f" AND protein_id IN (SELECT {column} FROM {table} WHERE {column} {in_range})"
    rows = connect(args).execute(
        "SELECT protein_id, sequence "
        f"FROM read_parquet('{args.release}/mgy_protein_sequences.parquet') WHERE {where}"
    )
    n, residues = 0, 0
    with open(args.out, "w") as out:
        while batch := rows.fetchmany(BATCH):
            n += len(batch)
            residues += sum(len(s) for _, s in batch)
            out.writelines(f">{i}\n{s[::-1] if args.reverse else s}\n" for i, s in batch)
    if n != args.proteins:
        sys.exit(f"[{args.lo}, {args.hi}): {n} sequences, planned {args.proteins}")
    Path(f"{args.out}.tsv").write_text(f"proteins\tresidues\n{n}\t{residues}\n")


def lines_of(path: str) -> Iterator[str]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        yield from f


def read_domtbl(path: str) -> pl.DataFrame:
    columns: dict[str, list[float]] = {name: [] for name in SCHEMA}
    for line in lines_of(path):
        if line.startswith("#"):
            continue
        c = line.split(maxsplit=22)
        c[4] = accession(c[4])
        for name, i in zip(SCHEMA, DOMTBL_COLUMNS, strict=True):
            columns[name].append(c[i])
    return pl.DataFrame(
        {k: pl.Series(k, v, strict=False).cast(SCHEMA[k]) for k, v in columns.items()},
        schema=SCHEMA,
    )


def cigar(model: str, target: str) -> tuple[str, int, int]:
    """One alignment's operations, and the model positions and residues they use."""
    ops: list[list] = []
    used = {"M": 0, "I": 0, "D": 0}
    for m, t in zip(model, target, strict=True):
        op = "I" if m == "." else "D" if t == "-" else "M"
        used[op] += 1
        if ops and ops[-1][0] == op:
            ops[-1][1] += 1
        else:
            ops.append([op, 1])
    return "".join(f"{n}{op}" for op, n in ops), used["M"] + used["D"], used["M"] + used["I"]


def read_alignments(path: str) -> pl.DataFrame:
    """Each domain's alignment from ``hmmsearch --notextw`` text output."""
    rows: list[tuple[int, int, int, str, str]] = []
    lines = lines_of(path)
    acc = None
    for line in lines:
        if line.startswith("Accession:"):
            acc = accession(line.split()[1])
        elif line.startswith("  == domain"):
            domain = int(line.split()[2])
            for line in lines:  # past any CS and RF lines, which end in a label, not a number
                model = line.split()
                if len(model) == 4 and model[1].isdigit() and model[3].isdigit():
                    break
            next(lines)  # the match line
            target = next(lines).split()
            pp = next(lines).split()[0]
            ops, positions, residues = cigar(model[2], target[2])
            if (positions, residues) != (
                int(model[3]) - int(model[1]) + 1,
                int(target[3]) - int(target[1]) + 1,
            ):
                sys.exit(f"{path}: {target[0]} domain {domain}: {ops} disagrees with coordinates")
            rows.append((int(target[0]), acc, domain, ops, pp))
    schema = {
        "protein_id": pl.Int64,
        "pfam_accession": pl.Int32,
        "domain": pl.Int16,
        "cigar": pl.String,
        "pp": pl.String,
    }
    return pl.DataFrame(rows, schema=schema, orient="row")


def parse(args: argparse.Namespace) -> None:
    thresholds = pl.read_csv(
        args.ga,
        separator="\t",
        schema_overrides={"pfam_accession": pl.Int32, "ga_seq": pl.Float32, "ga_dom": pl.Float32},
    )
    domains = read_domtbl(args.domtbl)
    keys = ["protein_id", "pfam_accession", "domain"]
    alignments = read_alignments(args.alignments)
    # Every domain has one alignment, and they number domains alike: a mismatch is a parse error.
    if domains.height != alignments.height:
        sys.exit(f"{args.domtbl}: {domains.height} domains, {alignments.height} alignments")
    hits = domains.join(alignments, on=keys, how="inner", validate="1:1")
    if hits.height != domains.height:
        sys.exit(f"{args.domtbl}: {domains.height - hits.height} domains with no alignment")
    if not args.pp:
        hits = hits.drop("pp")
    hits = hits.join(thresholds, on="pfam_accession", how="left")
    if hits["ga_seq"].null_count():
        sys.exit(f"{args.domtbl}: models with no GA threshold in {args.ga}")
    hits = hits.with_columns(
        above_ga=(pl.col("seq_score") >= pl.col("ga_seq")) & (pl.col("score") >= pl.col("ga_dom"))
    ).drop("ga_seq", "ga_dom")
    hits = hits.sort("protein_id", "ali_from")  # hmmsearch's order is by model
    hits.write_parquet(args.out, compression="zstd")
    counts = pl.read_csv(args.counts, separator="\t").row(0, named=True)
    above = hits.filter("above_ga")
    row = {
        "arm": args.arm,
        **counts,
        "domains": hits.height,
        "domains_above_ga": above.height,
        "proteins_with_hit": hits["protein_id"].n_unique(),
        "proteins_above_ga": above["protein_id"].n_unique(),
    }
    pl.DataFrame([row]).write_csv(args.stats, separator="\t")


def summary(args: argparse.Namespace) -> None:
    total = (
        pl.concat(pl.read_csv(p, separator="\t") for p in args.stats)
        .group_by("arm")
        .sum()
        .sort("arm")
    )
    total.write_csv(args.out, separator="\t")
    with open(args.out) as f:
        print(*csv.reader(f, delimiter="\t"), sep="\n")


def compare(args: argparse.Namespace) -> None:
    """Hit files are named ``{arm}.chunkN.parquet``."""
    tables = {
        arm: pl.scan_parquet([p for p in args.hits if Path(p).name.startswith(f"{arm}.")])
        .group_by((pl.col("score") / SCORE_BIN).floor().cast(pl.Int32) * SCORE_BIN, "above_ga")
        .len(arm)
        .collect()
        for arm in ("real", "reversed")
    }
    table = (
        tables["real"]
        .join(tables["reversed"], on=["score", "above_ga"], how="full", coalesce=True)
        .fill_null(0)
        .rename({"score": "score_from"})
        .sort("above_ga", "score_from")
    )
    table.write_csv(args.out, separator="\t")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--memory", default="4GB")
    sub = parser.add_subparsers(dest="step", required=True)
    p = sub.add_parser("ga")
    p.add_argument("--hmm", required=True)
    p.add_argument("--out", default="ga.tsv")
    for name in ("plan", "fasta"):
        p = sub.add_parser(name)
        p.add_argument("--release", required=True)
        p.add_argument("--targets", choices=["reps", "all"], default="reps")
        p.add_argument("--sample", type=int, default=1)
    p = sub.choices["plan"]
    p.add_argument("--chunk-proteins", type=int, default=500_000)
    p.add_argument("--bin-width", type=int, default=BIN_WIDTH)
    p.add_argument("--out", default="ranges.tsv")
    p = sub.choices["fasta"]
    p.add_argument("--lo", type=int, required=True)
    p.add_argument("--hi", type=int, required=True)
    p.add_argument("--proteins", type=int, required=True)
    p.add_argument("--reverse", action="store_true")
    p.add_argument("--out", required=True)
    p = sub.add_parser("parse")
    p.add_argument("--arm", default="real")
    p.add_argument("--domtbl", required=True)
    p.add_argument("--alignments", required=True, help="hmmsearch --notextw text output")
    p.add_argument("--pp", action="store_true", help="keep posterior probabilities")
    p.add_argument("--ga", required=True)
    p.add_argument("--counts", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--stats", required=True)
    p = sub.add_parser("summary")
    p.add_argument("stats", nargs="+")
    p.add_argument("--out", default="summary.tsv")
    p = sub.add_parser("compare")
    p.add_argument("hits", nargs="+")
    p.add_argument("--out", default="compare.tsv")
    args = parser.parse_args()
    steps = {
        "ga": ga,
        "plan": plan,
        "fasta": fasta,
        "parse": parse,
        "summary": summary,
        "compare": compare,
    }
    steps[args.step](args)


if __name__ == "__main__":
    main()
