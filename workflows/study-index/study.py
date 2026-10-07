"""Tables for a study index: members from MMseqs2 clusters, Pfam from hmmsearch."""

import argparse
import gzip
from collections.abc import Iterator
from pathlib import Path

import polars as pl


def fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """(header, sequence) records of a FASTA file, plain or gzip."""
    opener = gzip.open if str(path).endswith(".gz") else open
    header, chunks = None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def members(faa: list[str], clusters: str) -> pl.DataFrame:
    """The members table ``index`` reads. ``protein_id`` is the header's first word;
    ``full_length`` is False only for genes Prodigal marks partial (``partial=`` other than
    ``00``), so protein FASTA without that field counts as full length. ``clusters`` is
    ``mmseqs createtsv``'s (representative, member) pairs; a trailing ``*`` (Prodigal's
    stop) is dropped from sequences."""
    rows = [
        (h.split()[0], "partial=" not in h or "partial=00" in h, s.rstrip("*"))
        for path in faa
        for h, s in fasta(path)
    ]
    proteins = pl.DataFrame(
        rows, schema=["protein_id", "full_length", "sequence"], orient="row"
    ).unique("protein_id", keep="first")
    pairs = pl.read_csv(
        clusters, separator="\t", has_header=False, new_columns=["cluster_rep", "protein_id"]
    )
    table = pairs.join(proteins, on="protein_id", how="right").select(
        "protein_id",
        pl.col("cluster_rep").fill_null(pl.col("protein_id")),  # unclustered: singleton
        "full_length",
        "sequence",
    )
    return table.sort("cluster_rep", "protein_id")


def pfam(domtbl: list[str]) -> pl.DataFrame:
    """(``protein_id``, ``pfam_accession``) pairs from hmmsearch ``--domtblout`` files, the
    accession without its version."""
    rows = [
        (f[0], f[4].split(".")[0])
        for path in domtbl
        for line in Path(path).read_text().splitlines()
        if not line.startswith("#")
        for f in [line.split()]
    ]
    schema = {"protein_id": pl.String, "pfam_accession": pl.String}
    return pl.DataFrame(rows, schema=schema, orient="row").unique().sort(list(schema))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("members")
    p.add_argument("--faa", nargs="+", required=True)
    p.add_argument("--clusters", required=True)
    p.add_argument("--out", default="members.parquet")
    p = sub.add_parser("pfam")
    p.add_argument("--domtbl", nargs="+", required=True)
    p.add_argument("--out", default="pfam.parquet")
    args = parser.parse_args()
    if args.command == "members":
        members(args.faa, args.clusters).write_parquet(args.out)
    else:
        pfam(args.domtbl).write_parquet(args.out)


if __name__ == "__main__":
    main()
