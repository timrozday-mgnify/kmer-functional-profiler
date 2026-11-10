"""Read recruitment to confirmed targets, a step of its own (docs/tree-index-plan.md, Progress log:
recruitment on its own, 2026-11-10): the targets come from any tool that confirms them, such as
``compact.py query`` (the confirmed units' representatives) or kfpfam.

Every read pair is routed by seed votes on one diagonal band to its ``--route-top`` targets
(``_core.route_reads``, protein ``--route-k``-mers), joined with the ``--seeded`` pairs, then
pseudo-aligned (``_core.pseudo_align``'s anchored scan, ungapped BLOSUM62 bits). Each read pair
keeps its hits within ``--ambiguous-bits`` of its best, is called above a read bar from decoys
(the reads reversed, routed and aligned the same way; ``query.decoy_bar`` at ``--share`` and
``--chance-scale``), and its weight is shared among its hits by EM.

Inputs: ``--reads`` (R1 R2 FASTQ), ``--targets`` (FASTA: id and sequence), optionally
``--seeded`` (TSV: ``read``, the read pair's name as in the FASTQ, and ``target``, a target id;
pairs another tool already holds). Writes to ``--out``: ``reads.parquet`` (read: its index in
the FASTQ, name, target, score, weight), ``targets.tsv`` (target, reads: EM weight) and
``recruit.json``.
"""

import argparse
import gzip
import json
import time
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import polars as pl
import query
import tier2
import verify

from kmer_functional_profiler import _core

EMPTY = {"read": pl.Int64, "target": pl.Int64, "score": pl.Float64}


def batches(r1: str, r2: str, n: int) -> Iterator[tuple[int, list[bytes], list[bytes]]]:
    """(first read pair, mates 1, mates 2), ``n`` read pairs at a time."""
    lo, m1, m2 = 0, [], []
    with gzip.open(r1, "rb") as f, gzip.open(r2, "rb") as g:
        for i, (x, y) in enumerate(zip(f, g, strict=True)):
            if i % 4 == 1:
                m1.append(x.strip())
                m2.append(y.strip())
                if len(m1) == n:
                    yield lo, m1, m2
                    lo, m1, m2 = lo + n, [], []
    if m1:
        yield lo, m1, m2


def read_names(r1: str) -> list[str]:
    """Each read pair's name: R1's first word, less a trailing /1."""
    with gzip.open(r1, "rt") as f:
        names = [x[1:].split()[0] for i, x in enumerate(f) if i % 4 == 0]
    return [n[:-2] if n.endswith("/1") else n for n in names]


def aligned(
    pairs: pl.DataFrame, frames: dict[int, bytes], seqs: dict[int, bytes], batch: int
) -> pl.DataFrame:
    """(read, target, score: bits) by the anchored pseudo-alignment."""
    if not pairs.height:
        return pl.DataFrame(schema=EMPTY)
    got = tier2.align_chunk(
        pairs.select("read", pos="target"), frames, seqs, verify.blosum62(), batch, "anchored"
    )
    return got.rename({"pos": "target"})


def routed(
    reads: list[str], seqs: dict[int, bytes], seeded: pl.DataFrame, a: argparse.Namespace
) -> pl.DataFrame:
    """Every read pair's hits within ``--ambiguous-bits`` of its best: routed by seed votes to
    the targets, with its seeded pairs, pseudo-aligned."""
    targets = np.array(sorted(seqs), dtype=np.int64)
    rep_list = [seqs[t] for t in targets.tolist()]
    k = a.route_k or tier2.route_k(sum(map(len, rep_list)))
    votes = 4 if k <= 5 else 3
    out = [pl.DataFrame(schema=EMPTY)]
    for lo, m1, m2 in batches(*reads, a.read_batch):
        hi = lo + len(m1)
        fr = _core.pair_frames(m1, m2)
        r, leaf, _ = _core.route_reads(
            fr, rep_list, k, verify.protein_classes(), band=8, min_votes=votes, top=a.route_top,
            max_occ=2000,
        )  # fmt: skip
        pairs = pl.concat([
            pl.DataFrame({"read": r.astype(np.int64) + lo, "target": targets[leaf]}),
            seeded.filter(pl.col("read").is_between(lo, hi - 1)).select("read", "target"),
        ]).unique()  # fmt: skip
        got = aligned(pairs, dict(zip(range(lo, hi), fr, strict=True)), seqs, a.batch_pairs)
        top = pl.col("score").max().over("read")
        out.append(got.filter(pl.col("score") >= top - a.ambiguous_bits))
    return pl.concat(out)


def fasta(path: str) -> tuple[list[str], list[bytes]]:
    ids, seqs = [], []
    for rec in Path(path).read_text().split(">")[1:]:
        head, _, body = rec.partition("\n")
        ids.append(head.split()[0])
        seqs.append(body.replace("\n", "").encode())
    return ids, seqs


def main() -> None:
    p = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    p.add_argument("--reads", nargs=2, required=True)
    p.add_argument("--targets", required=True, help="FASTA of the confirmed targets")
    p.add_argument("--seeded", help="TSV: read (name), target (id): pairs already held")
    p.add_argument("--out", required=True)
    p.add_argument("--share", type=float, default=0.01, help="decoy share at the read bar")
    p.add_argument("--chance-scale", type=float, default=1.0, help="query.decoy_bar's scale")
    p.add_argument("--ambiguous-bits", type=float, default=3.0)
    p.add_argument("--route-k", type=int, default=0, help="0: from the targets' residues")
    p.add_argument("--route-top", type=int, default=16)
    p.add_argument("--read-batch", type=int, default=200_000)
    p.add_argument("--batch-pairs", type=int, default=2_000_000)
    p.add_argument("--em-rounds", type=int, default=20)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    ids, raw = fasta(a.targets)
    seqs = {i: verify.mask(s) for i, s in enumerate(raw)}
    names = read_names(a.reads[0])
    seeded = pl.DataFrame(schema={"read": pl.Int64, "target": pl.Int64})
    if a.seeded:  # names and ids to indices; pairs naming unknown reads or targets drop out
        by_name = pl.DataFrame({"read": names, "r": np.arange(len(names))})
        by_id = pl.DataFrame({"target": ids, "t": np.arange(len(ids))})
        s = pl.read_csv(a.seeded, separator="\t", schema={"read": pl.String, "target": pl.String})
        seeded = s.join(by_name, on="read").join(by_id, on="target").select(read="r", target="t")
    info: dict[str, float | int] = {"targets": len(ids), "read_pairs": len(names),
                                    "seeded_pairs": seeded.height}  # fmt: skip
    decoy = query.reversed_reads(a.reads, out)
    calls = {"sample": routed(a.reads, seqs, seeded, a) if seqs else pl.DataFrame(schema=EMPTY)}
    calls["decoy"] = routed(decoy, seqs, seeded.clear(), a) if seqs else pl.DataFrame(schema=EMPTY)
    per_read = {n: c.group_by("read").agg(pl.col("score").max()) for n, c in calls.items()}
    bar, _ = query.decoy_bar(per_read["sample"], per_read["decoy"], a.share, scale=a.chance_scale)
    called = calls["sample"].filter(pl.col("score").max().over("read") >= bar)
    w = tier2.em_weights(called.rename({"target": "pos"}), a.em_rounds).rename({"pos": "target"})
    w = w.with_columns(
        name=pl.Series(names).gather(w["read"]), target=pl.Series(ids).gather(w["target"])
    ).select("read", "name", "target", "score", "weight")
    w.sort("read", "target").write_parquet(out / "reads.parquet")
    (
        pl.DataFrame({"target": ids})
        .join(w.group_by("target").agg(reads=pl.col("weight").sum()), on="target", how="left")
        .with_columns(pl.col("reads").fill_null(0.0))
        .write_csv(out / "targets.tsv", separator="\t", float_precision=3)
    )
    for f in decoy:
        Path(f).unlink()
    info |= {"read_bar": bar, "called_reads": called["read"].n_unique(),
             "recruit_s": time.perf_counter() - t0, "peak_rss_gb": tier2.peak_rss_gb()}  # fmt: skip
    (out / "recruit.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info))


if __name__ == "__main__":
    main()
