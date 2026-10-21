"""What routing drops, searched against every profile unit (plan, phase 13, routing).

``dropped``: the pairs routing (k = 8, min_kmers = 3 by default) does not send to DIAMOND,
as FASTA (names cut at the first word), a seeded subsample of ``--pairs`` of them
(0 = all), and ``dropped.json`` (pairs, dropped, sampled).

``fullref``: one FASTA of every unit in the profile, ``--members`` sequences each, and a
table of each sequence's unit, its role in the per-sample reference (target, neighbour,
``profile``: in the profile only, or ``withheld``: left out on purpose, the spike's
siblings) and its ``present_q``.

``summarise``: per arm (a full reference; DIAMOND tabular ``qseqid sseqid pident length
bitscore`` per mate), each dropped pair's best alignment over both mates: pairs aligned, by
identity bin and by the role of the unit hit; high-identity pairs (>= ``HIGH``) on targets
or neighbours (routing's misses) and on profile-only units (reads of units the reference
lacks), and the profile-only units they gather on; with ``--assignments``, how the unrouted
baseline treated the same pairs. Between two arms, the identity each pair gains with more
members. With ``--trace``, DIAMOND's CPU-h per arm and that per million dropped pairs.
"""

import argparse
import gzip
import json
import random
import tempfile
from pathlib import Path

import numpy as np
import polars as pl

from kfp_pileup.reference import REF_FASTA, REFS
from kfp_pileup.route import MIN_KMERS, ROUTE_K, _records, blocks, kmer_table, pair_hits
from kmer_functional_profiler.index import read_fasta
from kmer_functional_profiler.query import fdr_q
from kmer_functional_profiler.store import fetch

HIGH = 90.0  # percent identity: a read of the unit itself, or a close strain
BINS = [50.0, 70.0, 90.0, 95.0]
FIELDS = ["qseqid", "sseqid", "pident", "length", "bitscore"]
MATE = r"/[12]$"


def dropped(args: argparse.Namespace) -> None:
    params = json.loads((args.index / "meta.json").read_text())["params"]
    refs = pl.read_csv(args.ref / REFS, separator="\t", infer_schema_length=None)
    sequences = dict(read_fasta(args.ref / REF_FASTA))
    seqs = [sequences[s].encode() for s in refs["seq_id"]]
    block = blocks(refs, 2_000_000)
    args.out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.out) as scratch:
        table = kmer_table(seqs, block, args.k, params["alphabet"], Path(scratch))
        rows, pairs = pair_hits(table, args.reads, args.k, params["alphabet"])
        del table
    routed = set(rows.filter(pl.col("kmers") >= args.min_kmers)["read"].unique().to_list())
    drop = [i for i in range(pairs) if i not in routed]
    sample = (
        sorted(random.Random(1).sample(drop, args.pairs)) if 0 < args.pairs < len(drop) else drop
    )
    keep = np.zeros(pairs, dtype=bool)
    keep[sample] = True
    for mate, path in enumerate(args.reads, 1):
        with gzip.open(args.out / f"dropped_{mate}.fa.gz", "wb", compresslevel=1) as f:
            for i, (name, seq) in enumerate(_records(path)):
                if keep[i]:
                    f.write(b">" + name + b"\n" + seq + b"\n")
    counts = {"pairs": pairs, "routed": len(routed), "dropped": len(drop), "sampled": len(sample)}
    counts |= {"k": args.k, "min_kmers": args.min_kmers}
    (args.out / "dropped.json").write_text(json.dumps(counts, indent=1) + "\n")
    print(counts)


def fullref(args: argparse.Namespace) -> None:
    with open(args.profile) as f:
        has_q = "present_q" in f.readline().rstrip("\n").split("\t")
    prof = pl.read_csv(
        args.profile,
        separator="\t",
        columns=["unit", "present_prob"] + (["present_q"] if has_q else []),
    ).with_columns(pl.col("unit").cast(pl.UInt32))
    if not has_q:
        prof = prof.with_columns(present_q=pl.Series(fdr_q(prof["present_prob"].to_numpy())))
    roles = (
        pl.read_csv(args.ref / REFS, separator="\t", columns=["unit", "role"])
        .with_columns(pl.col("unit").cast(pl.UInt32))
        .unique("unit")
    )
    held = args.withhold.read_text().split() if args.withhold else []
    seqs = fetch(args.index, prof["unit"].to_list(), args.members, args.store or ())
    seqs = (
        seqs.with_columns(pl.col("unit").cast(pl.UInt32))
        .join(roles, on="unit", how="left")
        .join(prof.select("unit", "present_q"), on="unit", how="left")
        .with_columns(
            pl.when(pl.col("cluster_rep").cast(pl.String).is_in(held))
            .then(pl.lit("withheld"))
            .otherwise(pl.col("role").fill_null("profile"))
            .alias("role"),
            sseqid=pl.format("f{}", pl.int_range(pl.len())),
            length=pl.col("sequence").str.len_bytes(),
        )
    )
    with open(args.out_fasta, "w") as f:
        for sid, seq in seqs.select("sseqid", "sequence").iter_rows():
            f.write(f">{sid}\n{seq}\n")
    seqs.drop("sequence").write_csv(args.out_table, separator="\t")
    print({"units": prof.height, "sequences": seqs.height, "residues": int(seqs["length"].sum())})


def best_per_pair(r1: Path, r2: Path, table: Path) -> pl.DataFrame:
    """Each dropped pair's best alignment over both mates, with the unit hit and its role."""
    cols = pl.read_csv(table, separator="\t", columns=["sseqid", "unit", "role", "present_q"])
    hits = pl.concat(
        pl.read_csv(p, separator="\t", has_header=False, new_columns=FIELDS)
        for p in (r1, r2)
        if p.stat().st_size
    ) if any(p.stat().st_size for p in (r1, r2)) else pl.DataFrame(
        schema={"qseqid": pl.String, "sseqid": pl.String, "pident": pl.Float64,
                "length": pl.Int64, "bitscore": pl.Float64}
    )  # fmt: skip
    return (
        hits.with_columns(read_name=pl.col("qseqid").str.replace(MATE, ""))
        .sort("bitscore", descending=True)
        .group_by("read_name")
        .first()
        .join(cols, on="sseqid", how="left")
    )


def arm_summary(best: pl.DataFrame, n_dropped: int) -> dict[str, float | int]:
    high = best.filter(pl.col("pident") >= HIGH)
    out: dict[str, float | int] = {
        "dropped_pairs": n_dropped,
        "aligned": best.height / n_dropped if n_dropped else 0.0,
    }
    edges = [0.0, *BINS, 101.0]
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        part = best.filter((pl.col("pident") >= lo) & (pl.col("pident") < hi))
        out[f"identity_{lo:g}_{min(hi, 100):g}"] = part.height / n_dropped if n_dropped else 0.0
    for role in ("target", "neighbour", "profile", "withheld"):
        out[f"best_on_{role}"] = best.filter(pl.col("role") == role).height / max(n_dropped, 1)
        out[f"high_on_{role}"] = high.filter(pl.col("role") == role).height / max(n_dropped, 1)
    units = (
        high.filter(pl.col("role") == "profile")
        .group_by("unit")
        .agg(pairs=pl.len(), present_q=pl.col("present_q").first())
    )
    out["profile_units_high"] = units.height
    out["profile_units_high_10"] = units.filter(pl.col("pairs") >= 10).height
    out["profile_units_high_median_q"] = (
        float(units["present_q"].median()) if units.height else float("nan")  # type: ignore[arg-type]
    )
    return out


def summarise(args: argparse.Namespace) -> None:
    n = json.loads((args.dropped / "dropped.json").read_text())["sampled"]
    names = pl.Series(
        "read_name", [nm.decode() for nm, _ in _records(args.dropped / "dropped_1.fa.gz")]
    ).str.replace(MATE, "")
    rows, bests = [], {}
    for arm in args.arm:
        name, directory = arm.split("=", 1)
        d = Path(directory)
        bests[name] = best_per_pair(d / "r1.tsv", d / "r2.tsv", d / "full.tsv")
        rows.append({"arm": name, **arm_summary(bests[name], n)})
    table = pl.DataFrame(rows)
    if args.assignments:
        a = pl.read_parquet(args.assignments).with_columns(
            pl.col("read_name").str.replace(MATE, "")
        )
        dropped_a = a.join(pl.DataFrame(names), on="read_name", how="semi")
        assigned = (
            dropped_a.group_by("read_name", "mate", "unit")
            .agg(pl.col("em_posterior").sum())
            .filter(pl.col("em_posterior") >= 0.5)["read_name"]
            .n_unique()
        )
        novel = dropped_a.group_by("read_name").agg(pl.col("em_novel").max())
        table = table.with_columns(
            baseline_aligned=pl.lit(dropped_a["read_name"].n_unique() / max(n, 1)),
            baseline_unit=pl.lit(assigned / max(n, 1)),
            baseline_novel=pl.lit(novel.filter(pl.col("em_novel") >= 0.5).height / max(n, 1)),
        )
    if args.trace and args.trace.exists():
        t = pl.read_csv(args.trace, separator="\t", infer_schema_length=None, null_values=["-"])
        t = t.filter(pl.col("status").is_in(["COMPLETED", "CACHED"])).with_columns(
            cpu_h=pl.col("realtime").cast(pl.Float64)
            / 3.6e6
            * pl.col("%cpu").cast(pl.Float64)
            / 100,
            arm=pl.col("name").str.extract(r"DIAMOND \((\w+)\)"),
        )
        cost = (
            t.filter(pl.col("arm").is_not_null())
            .group_by("arm")
            .agg(
                diamond_cpu_h=pl.col("cpu_h").sum(),
                diamond_peak_rss_gb=pl.col("peak_rss").cast(pl.Float64).max() / 1e9,
            )
        )
        table = table.join(cost, on="arm", how="left").with_columns(
            cpu_h_per_m_pairs=pl.col("diamond_cpu_h") / max(n, 1) * 1e6
        )
    if len(bests) == 2:  # identity gained with more members, per pair aligned in either
        (a1, b1), (a2, b2) = bests.items()
        both = (
            b1.select("read_name", p1="pident")
            .join(b2.select("read_name", p2="pident"), on="read_name", how="full", coalesce=True)
            .fill_null(0.0)
        )
        gain = both["p2"] - both["p1"]
        print(
            f"{a2} vs {a1}: pairs aligned {b2.height} vs {b1.height}; identity gain median "
            f"{gain.median():.2f}, >= 2 points in {(gain >= 2).mean():.3f} of pairs aligned"
            f" in either; reaching >= {HIGH:g}% only with {a2}: "
            f"{((both['p2'] >= HIGH) & (both['p1'] < HIGH)).sum()} pairs"
        )
    table.write_csv(args.out, separator="\t")
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=250, fmt_float="mixed"):
        print(table.transpose(include_header=True, column_names="arm"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Routing's dropped pairs against every unit")
    sub = parser.add_subparsers(dest="command", required=True)
    d = sub.add_parser("dropped")
    d.add_argument("--ref", type=Path, required=True)
    d.add_argument("--index", type=Path, required=True)
    d.add_argument("--reads", type=Path, nargs=2, required=True)
    d.add_argument("--k", type=int, default=ROUTE_K)
    d.add_argument("--min-kmers", type=int, default=MIN_KMERS)
    d.add_argument("--pairs", type=int, default=1_000_000, help="subsample (0 = all)")
    d.add_argument("--out", type=Path, default=Path("."))
    f = sub.add_parser("fullref")
    f.add_argument("--profile", type=Path, required=True)
    f.add_argument("--ref", type=Path, required=True)
    f.add_argument("--index", type=Path, required=True)
    f.add_argument("--store", type=Path, nargs="*")
    f.add_argument("--members", type=int, default=1)
    f.add_argument("--withhold", type=Path, help="units left out of the reference (spike)")
    f.add_argument("--out-fasta", type=Path, default=Path("full.fa"))
    f.add_argument("--out-table", type=Path, default=Path("full.tsv"))
    s = sub.add_parser("summarise")
    s.add_argument("--dropped", type=Path, required=True, help="directory of `dropped`")
    s.add_argument("--arm", action="append", required=True, help="name=dir (r1/r2.tsv, full.tsv)")
    s.add_argument("--assignments", type=Path)
    s.add_argument("--trace", type=Path)
    s.add_argument("--out", type=Path, default=Path("summary.tsv"))
    args = parser.parse_args()
    {"dropped": dropped, "fullref": fullref, "summarise": summarise}[args.command](args)


if __name__ == "__main__":
    main()
