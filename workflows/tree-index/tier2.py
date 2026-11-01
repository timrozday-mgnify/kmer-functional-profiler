"""Tree index, tier 2: a per-sample database over tier 1's candidates, for read classification
(docs/tree-index-plan.md, Uses beyond tier 1; Progress log: tier 2 construction).

Tier 1's lessons carried over: an alignment, not k-mer counts, makes the call (``verify.py``'s
anchored ungapped scan, tandem repeats masked); group representatives keep recall, so tier 1
routes each read to groups; decoys (the reads reversed) set the bar; reads cannot tell
near-identical leaves apart, so ties are kept, not broken.

1. **Work list.** Tier 1's aligned (read, leaf) pairs (``verify.py``'s ``aligned.parquet``,
   ``aligned_decoy.parquet``) within ``--route-below`` bits of the read-level bar, each leaf
   expanded to every leaf of its layer group (``layer_group.npy``): a read reaches the member
   closest to it, not only its group's composite.
2. **Chunks along the tree.** Candidate leaves in leaf order (tier 1's merge-tree order, so
   relatives are adjacent and a chunk is a contiguous range) are packed into chunks of at most
   ``--chunk-residues`` residues. One streaming pass over ``--members`` writes the candidates'
   sequences (each leaf's longest protein) in leaf order (``seqs.parquet``), so a chunk reads
   one range of it. A layer group is a single-linkage cluster of the tree, so one contiguous
   range: the work list splits by chunk once, each group in one chunk or two.
3. **Alignment per chunk,** in batches of at most ``--batch-pairs`` pairs, each batch's scores
   spilled to ``--out``/spill/; then each read's best leaves (ties kept) and a read-level bar
   at each ``--shares`` where decoy reads called are at most that share of the sample's.

Memory: one chunk's sequences, one batch of pairs, and the routed reads' frames (~300 B per
read pair; read from ``verify.py``'s ``frames_*.parquet``, recomputed only if absent, which
costs every read's frames). The work list is counted before it is expanded, and expanded
chunk by chunk.

Writes to ``--out``: ``reads.parquet`` (read, pos, score: each read pair's best leaves),
``cand_verify_<share>.npy`` (leaves best for a read at the bar: tier 2's presence),
``verify.json`` (``read_bar_<share>``, as ``verify.py`` names it, so ``recall.py reads`` scores
tier 2 as it scores tier 1) and ``tier2.json`` (work list, chunks, times, peak RSS).
"""

import argparse
import json
import math
import resource
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
import query
import verify
from numpy.typing import NDArray

from kmer_functional_profiler import _core


def peak_rss_gb() -> float:
    """This process's peak RSS (macOS reports bytes, Linux KiB)."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / (2**30 if sys.platform == "darwin" else 2**20)


def candidates(index: Path, qdir: Path, share: str) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Tier 1's candidate leaves (``cand_verify_<share>``, each representative's group
    expanded) and every position's group."""
    pos_len = np.load(index / "pos_len.npy")
    gfile = index / "layer_group.npy"
    group = np.load(gfile).astype(np.int64) if gfile.exists() else np.arange(len(pos_len))
    cand = np.load(qdir / f"cand_verify_{share}.npy").astype(np.int64)
    return np.flatnonzero(np.isin(group, group[cand])), group


def chunks(leaves: NDArray[np.int64], length: NDArray[np.int64], budget: int) -> NDArray[np.int64]:
    """Each leaf's chunk: leaves in leaf order packed into runs of at most ``budget``
    residues (a leaf longer than that is a chunk of its own)."""
    out = np.zeros(len(leaves), dtype=np.int64)
    c, used = 0, 0
    for i, n in enumerate(length[leaves].tolist()):
        if used and used + n > budget:
            c, used = c + 1, 0
        out[i], used = c, used + n
    return out


def leaf_store(index: Path, members: str, leaves: NDArray[np.int64], dst: Path) -> None:
    """Each candidate leaf's longest protein (ties: the first) to ``dst``, sorted by
    position, in one streaming pass over the members (row groups prune a chunk's reads)."""
    prot = (
        pl.scan_parquet(index / "proteins.parquet")
        .filter(pl.col("pos").is_in(pl.Series(leaves).implode()))
        .select("protein_id", "pos")
    )
    (
        pl.scan_parquet(members)
        .select("protein_id", "sequence")
        .join(prot, on="protein_id")
        .with_columns(n=pl.col("sequence").str.len_bytes())
        .sort("pos", "n", descending=[False, True])
        .unique("pos", keep="first", maintain_order=True)
        .select("pos", "sequence")
        .sink_parquet(dst, row_group_size=20_000)
    )


def chunk_seqs(store: Path, lo: int, hi: int) -> dict[int, bytes]:
    """Leaves ``lo..hi``'s sequences from the leaf-ordered store, masked."""
    got = pl.scan_parquet(store).filter(pl.col("pos").is_between(lo, hi)).collect()
    return dict(
        zip(got["pos"].to_list(), (verify.mask(s.encode()) for s in got["sequence"]), strict=True)
    )


def work_list(aligned: pl.DataFrame, group: NDArray[np.int64], floor: float) -> pl.DataFrame:
    """(read, g): each read's groups with a tier-1 pair scoring at least ``floor``."""
    kept = aligned.filter(pl.col("score") >= floor)
    return kept.select("read", g=pl.Series(group[kept["pos"].to_numpy()])).unique()


def routed_frames(
    tier1: query.Index, saved: Path, reads: list[str], wl: pl.DataFrame
) -> dict[int, bytes]:
    """The routed reads' frames: from ``verify.py``'s ``frames_*.parquet`` (only those rows
    read), else recomputed (``verify.read_frames``, the memory peak: every read's frames)."""
    routed = wl["read"].unique()
    if saved.exists():
        got = pl.scan_parquet(saved).filter(pl.col("read").is_in(routed.implode())).collect()
        return dict(zip(got["read"].to_list(), got["frames"].to_list(), strict=True))
    frames = verify.read_frames(tier1, *reads)
    return {r: frames[r] for r in routed.to_list()}


def align_chunk(
    pairs: pl.DataFrame,
    frames: dict[int, bytes],
    seqs: dict[int, bytes],
    matrix: NDArray[np.int8],
    batch: int,
) -> pl.DataFrame:
    """(read, pos, score) for a chunk's pairs, ``batch`` pairs at a time: each batch passes
    only its own reads' frames and leaves' sequences to the kernel."""
    out = []
    pairs = pairs.sort("pos", "read")  # one sketch per run of a leaf
    for start in range(0, pairs.height, batch):
        b = pairs.slice(start, batch)
        reads, ri = np.unique(b["read"].to_numpy(), return_inverse=True)
        pos, pi = np.unique(b["pos"].to_numpy(), return_inverse=True)
        s = _core.pseudo_align(
            [frames[r] for r in reads], [seqs[p] for p in pos.tolist()],
            ri.astype(np.uint32), pi.astype(np.uint32), "anchored", 4, verify.protein_classes(),
            matrix, sketch=(8, 1, 8), band=8, top=2,
        )  # fmt: skip
        bits = (verify.LAMBDA * s.astype(np.float64) - math.log(verify.K)) / math.log(2)
        out.append(b.select("read", "pos").with_columns(score=bits))
    return pl.concat(out) if out else pairs.select("read", "pos").with_columns(score=pl.lit(0.0))


def main() -> None:
    p = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    p.add_argument("--index", required=True)
    p.add_argument("--members", required=True, help="the index's members (sequences)")
    p.add_argument("--reads", nargs=2, required=True)
    p.add_argument("--query", required=True, help="query.py and verify.py's --out")
    p.add_argument("--share", default="0.01", help="tier 1's candidate and read bar share")
    p.add_argument("--route-below", type=float, default=10.0,
                   help="route tier-1 pairs within this many bits of its read bar")  # fmt: skip
    p.add_argument("--chunk-residues", type=int, default=50_000_000)
    p.add_argument("--batch-pairs", type=int, default=2_000_000)
    p.add_argument("--shares", default="0.01,0.05,0.1")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    t0 = time.perf_counter()
    index, qdir, out = Path(a.index), Path(a.query), Path(a.out)
    spill = out / "spill"
    spill.mkdir(parents=True, exist_ok=True)
    tier1 = query.Index(a.index)
    matrix = verify.blosum62()
    leaves, group = candidates(index, qdir, a.share)
    length = np.load(index / "pos_len.npy").astype(np.int64)
    chunk_of = chunks(leaves, length, a.chunk_residues)
    leaf_chunk = pl.DataFrame({"pos": leaves, "g": group[leaves], "chunk": chunk_of})
    floor = json.loads((qdir / "verify.json").read_text())[f"read_bar_{a.share}"] - a.route_below
    info: dict[str, float | int] = {
        "candidate_leaves": len(leaves),
        "chunks": int(chunk_of.max()) + 1 if len(leaves) else 0,
        "residues": int(length[leaves].sum()),
        "route_floor": floor,
    }
    store = out / "seqs.parquet"
    leaf_store(index, a.members, leaves, store)
    info["store_s"] = time.perf_counter() - t0
    # chunks holding each group (one, or two where a chunk boundary cuts it)
    g_chunk = leaf_chunk.select("g", "chunk").unique()
    saved = (qdir / "frames_decoy.parquet").exists()  # else the decoy reads to recompute from
    reads = {"sample": a.reads, "decoy": [] if saved else query.reversed_reads(a.reads, out)}
    for name, aligned in (("sample", "aligned.parquet"), ("decoy", "aligned_decoy.parquet")):
        wl = work_list(pl.read_parquet(qdir / aligned), group, floor)
        # counted before expanding: (read, group) x the group's candidate leaves
        size = wl.join(leaf_chunk.group_by("g").agg(n=pl.len()), on="g")["n"].sum()
        info |= {f"{name}_reads_routed": wl["read"].n_unique(), f"{name}_pairs": int(size)}
        t1 = time.perf_counter()
        frames = routed_frames(tier1, qdir / f"frames_{name}.parquet", reads[name], wl)
        info[f"{name}_frames_s"] = time.perf_counter() - t1
        t1 = time.perf_counter()
        parts = wl.join(g_chunk, on="g").partition_by("chunk", as_dict=True)
        for (c,), part in parts.items():
            mine = leaf_chunk.filter(pl.col("chunk") == c)
            pairs = part.join(mine.select("g", "pos"), on="g").select("read", "pos")
            seqs = chunk_seqs(store, int(mine["pos"].min()), int(mine["pos"].max()))
            got = align_chunk(pairs, frames, seqs, matrix, a.batch_pairs)
            got.write_parquet(spill / f"{name}.{c:05d}.parquet")
        info[f"{name}_align_s"] = time.perf_counter() - t1
        del frames
    # each read's best leaves, merged across chunks lazily
    best = {}
    schema = {"read": pl.Int64, "pos": pl.Int64, "score": pl.Float64}
    for name in ("sample", "decoy"):
        parts = sorted(spill.glob(f"{name}.*.parquet"))
        if not parts:  # nothing routed (a small or clean sample's decoys)
            best[name] = pl.DataFrame(schema=schema)
            continue
        lf = pl.scan_parquet(parts).cast(schema)  # type: ignore[arg-type]
        best[name] = lf.filter(pl.col("score") == pl.col("score").max().over("read")).collect()
    best["sample"].sort("read", "pos").write_parquet(out / "reads.parquet")
    per_read = {n: b.group_by("read").agg(pl.col("score").first()) for n, b in best.items()}
    bars = {}
    for s in a.shares.split(","):
        bar = query.decoy_bar(per_read["sample"], per_read["decoy"], float(s))[0]
        bars[f"read_bar_{float(s)}"] = bar
        called = best["sample"].filter(pl.col("score") >= bar)["pos"].unique().sort()
        np.save(out / f"cand_verify_{float(s)}.npy", called.to_numpy())
    (out / "verify.json").write_text(json.dumps(bars, indent=1) + "\n")
    shutil.rmtree(spill)
    store.unlink()
    for f in out.glob("decoy_*.fq.gz"):
        f.unlink()
    info |= bars | {"total_s": time.perf_counter() - t0, "peak_rss_gb": peak_rss_gb()}
    (out / "tier2.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info))


if __name__ == "__main__":
    main()
