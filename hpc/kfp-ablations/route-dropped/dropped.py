"""What routing drops, searched against every profile unit (plan, phase 13, routing).

Reads are single-end here (the user's choice, 2026-10-21): each mate is routed on its own
and dropped when it shares fewer than ``min_kmers`` distinct k-mers with every block.

``dropped``: the dropped reads (k = 8, min_kmers = 3 by default) as one FASTA (``<pair
name>/<mate>``), a seeded subsample of ``--reads`` of them (0 = all); ``dropped.parquet``
(each sampled read's ``total_kmers``, its distinct k-mers in all, and ``shared_kmers``, the
most it shares with one block); ``dropped.json``: reads, routed, dropped (also with the
partner routed, which the pair-routed pipeline still aligns), and the dropped reads by
total k-mers (few: a short, N-rich or, with a quality mask, low-quality read; many but few
shared: a read the reference lacks or holds distantly).

``fullref``: one FASTA of every unit in the profile, ``--members`` sequences each, and a
table of each sequence's unit, its role in the per-sample reference (target, neighbour,
``profile``: in the profile only, or ``withheld``: left out on purpose, the spike's
siblings) and its ``present_q``.

``summarise``: per arm (a full reference; DIAMOND tabular ``qseqid sseqid pident length
bitscore``), each dropped read's best alignment: reads aligned, by identity bin and by the
role of the unit hit; high-identity reads (>= ``HIGH``) on targets or neighbours (routing's
misses) and on profile-only units (reads of units the reference lacks), and the units they
gather on; the same by total k-mers; with ``--assignments``, how the unrouted baseline
treated the same reads. Between two arms, the identity each read gains with more members.
With ``--trace``, DIAMOND's CPU-h per arm and that per million dropped reads.
"""

import argparse
import gzip
import itertools
import json
import re
import tempfile
from pathlib import Path

import numpy as np
import polars as pl

from kfp_pileup.reference import REF_FASTA, REFS
from kfp_pileup.route import (
    MIN_KMERS,
    ROUTE_K,
    TOTAL_BINS,
    _records,
    blocks,
    by_total,
    kmer_table,
    pair_hits,
)
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
        rows, totals, n = pair_hits(table, args.reads, args.k, params["alphabet"], single=True)
        del table
    total, shared = np.zeros(n, dtype=np.int64), np.zeros(n, dtype=np.int64)
    total[totals["read"].to_numpy().astype(np.int64)] = totals["total_kmers"].to_numpy()
    most = rows.group_by("read").agg(pl.col("kmers").max())
    shared[most["read"].to_numpy().astype(np.int64)] = most["kmers"].to_numpy()
    routed = shared >= args.min_kmers
    drop = np.flatnonzero(~routed)
    sample = drop
    if 0 < args.reads_n < len(drop):
        sample = np.sort(np.random.default_rng(1).choice(drop, args.reads_n, replace=False))
    keep = np.zeros(n, dtype=bool)
    keep[sample] = True
    with gzip.open(args.out / "dropped.fa.gz", "wb", compresslevel=1) as f:
        for mate, path in enumerate(args.reads):
            for i, (name, seq) in enumerate(_records(path)):
                if keep[2 * i + mate]:
                    stem = re.sub(rb"/[12]$", b"", name)
                    f.write(b">" + stem + b"/%d\n" % (mate + 1) + seq + b"\n")
    pl.DataFrame(
        {
            "read": sample,
            "mate": (sample % 2 + 1).astype(np.uint8),
            "total_kmers": total[sample],
            "shared_kmers": shared[sample],
        }
    ).write_parquet(args.out / "dropped.parquet")
    counts: dict[str, object] = {
        "reads": n,
        "routed": int(routed.sum()),
        "dropped": len(drop),
        "dropped_partner_routed": int(routed[drop ^ 1].sum()),
        "pairs_both_dropped": int((~routed[0::2] & ~routed[1::2]).sum()),
        "sampled": len(sample),
        "k": args.k,
        "min_kmers": args.min_kmers,
        "dropped_by_total_kmers": by_total(pl.Series(total[drop]), 0),
        "total_kmers_median": {
            "routed": float(np.median(total[routed])) if routed.any() else None,
            "dropped": float(np.median(total[drop])) if len(drop) else None,
        },
    }
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


def best_per_read(hits_path: Path, table: Path) -> pl.DataFrame:
    """Each dropped read's best alignment, with the unit hit and its role."""
    cols = pl.read_csv(table, separator="\t", columns=["sseqid", "unit", "role", "present_q"])
    if not hits_path.stat().st_size:
        return pl.DataFrame(
            schema={"qseqid": pl.String, "sseqid": pl.String, "pident": pl.Float64,
                    "length": pl.Int64, "bitscore": pl.Float64, "unit": pl.Int64,
                    "role": pl.String, "present_q": pl.Float64}
        )  # fmt: skip
    hits = pl.read_csv(hits_path, separator="\t", has_header=False, new_columns=FIELDS)
    return (
        hits.sort("bitscore", descending=True)
        .group_by("qseqid")
        .first()
        .join(cols, on="sseqid", how="left")
    )


def arm_summary(best: pl.DataFrame, n_dropped: int) -> dict[str, float | int]:
    high = best.filter(pl.col("pident") >= HIGH)
    out: dict[str, float | int] = {
        "dropped_reads": n_dropped,
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
        .agg(reads=pl.len(), present_q=pl.col("present_q").first())
    )
    out["profile_units_high"] = units.height
    out["profile_units_high_10"] = units.filter(pl.col("reads") >= 10).height
    out["profile_units_high_median_q"] = (
        float(units["present_q"].median()) if units.height else float("nan")  # type: ignore[arg-type]
    )
    return out


def summarise(args: argparse.Namespace) -> None:
    n = json.loads((args.dropped / "dropped.json").read_text())["sampled"]
    with gzip.open(args.dropped / "dropped.fa.gz", "rt") as f:  # in dropped.parquet's order
        names = pl.Series("qseqid", [line[1:].split()[0] for line in f if line.startswith(">")])
    info = pl.read_parquet(args.dropped / "dropped.parquet").with_columns(qseqid=names)
    edges = [0, *TOTAL_BINS, None]
    rows, bests = [], {}
    for arm in args.arm:
        name, directory = arm.split("=", 1)
        d = Path(directory)
        best = best_per_read(d / "r.tsv", d / "full.tsv")
        bests[name] = best
        rows.append({"arm": name, "total_kmers": "all", **arm_summary(best, n)})
        for lo, hi in itertools.pairwise(edges):  # the same by the reads' total k-mers
            t = pl.col("total_kmers")
            part = info.filter(t >= lo if hi is None else (t >= lo) & (t < hi))
            label = f"{lo}+" if hi is None else f"{lo}-{hi - 1}"
            sub = best.join(part.select("qseqid"), on="qseqid", how="semi")
            rows.append({"arm": name, "total_kmers": label, **arm_summary(sub, part.height)})
    table = pl.DataFrame(rows)
    if args.assignments:
        a = pl.read_parquet(args.assignments).with_columns(
            qseqid=pl.format("{}/{}", pl.col("read_name").str.replace(MATE, ""), pl.col("mate") + 1)
        )
        dropped_a = a.join(pl.DataFrame(names), on="qseqid", how="semi")
        assigned = (
            dropped_a.group_by("qseqid", "unit")
            .agg(pl.col("em_posterior").sum())
            .filter(pl.col("em_posterior") >= 0.5)["qseqid"]
            .n_unique()
        )
        novel = dropped_a.group_by("qseqid").agg(pl.col("em_novel").max())
        table = table.with_columns(
            baseline_aligned=pl.lit(dropped_a["qseqid"].n_unique() / max(n, 1)),
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
            cpu_h_per_m_reads=pl.col("diamond_cpu_h") / max(n, 1) * 1e6
        )
    if len(bests) == 2:  # identity gained with more members, per read aligned in either
        (a1, b1), (a2, b2) = bests.items()
        both = (
            b1.select("qseqid", p1="pident")
            .join(b2.select("qseqid", p2="pident"), on="qseqid", how="full", coalesce=True)
            .fill_null(0.0)
        )
        gain = both["p2"] - both["p1"]
        print(
            f"{a2} vs {a1}: reads aligned {b2.height} vs {b1.height}; identity gain median "
            f"{gain.median():.2f}, >= 2 points in {(gain >= 2).mean():.3f} of reads aligned"
            f" in either; reaching >= {HIGH:g}% only with {a2}: "
            f"{((both['p2'] >= HIGH) & (both['p1'] < HIGH)).sum()} reads"
        )
    table.write_csv(args.out, separator="\t")
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=250, fmt_float="mixed"):
        print(table)


def main() -> None:
    parser = argparse.ArgumentParser(description="Routing's dropped pairs against every unit")
    sub = parser.add_subparsers(dest="command", required=True)
    d = sub.add_parser("dropped")
    d.add_argument("--ref", type=Path, required=True)
    d.add_argument("--index", type=Path, required=True)
    d.add_argument("--reads", type=Path, nargs=2, required=True)
    d.add_argument("--k", type=int, default=ROUTE_K)
    d.add_argument("--min-kmers", type=int, default=MIN_KMERS)
    d.add_argument("--reads-n", type=int, default=2_000_000, help="dropped reads aligned (0 = all)")
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
    s.add_argument("--arm", action="append", required=True, help="name=dir (r.tsv, full.tsv)")
    s.add_argument("--assignments", type=Path)
    s.add_argument("--trace", type=Path)
    s.add_argument("--out", type=Path, default=Path("summary.tsv"))
    args = parser.parse_args()
    {"dropped": dropped, "fullref": fullref, "summarise": summarise}[args.command](args)


if __name__ == "__main__":
    main()
