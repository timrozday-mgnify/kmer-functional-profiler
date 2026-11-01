"""Tree index, tier 2: a per-sample database over tier 1's candidates, for read classification
(docs/tree-index-plan.md, Uses beyond tier 1; Progress log: tier 2 construction).

Tier 1's lessons carried over: an alignment, not k-mer counts, makes the call (``verify.py``'s
anchored scan, tandem repeats masked); group representatives keep recall, so tier 1 routes
each read to groups; decoys (the reads reversed) set the bar; reads cannot tell near-identical
leaves apart, so ties are shared, not broken.

1. **Candidates and routing.** The candidates are every leaf of tier 1's present groups
   (``cand_verify_<share>``) and of the groups its aligned pairs route reads to (within
   ``--route-below`` bits of its read bar), with, at ``--neighbours`` d, the groups next to
   them on the spanning tree (an edge at dense-set distance at most d). Tier 1 is a profile:
   its sparse keys seed only some of a strain's reads, so ``--route dense`` (the default)
   routes every read pair itself: all its twelve frames vote, by protein ``--route-k``-mer
   seeds on one diagonal band (``_core.route_reads``, the candidates' k-mers held more than
   ``--route-max-occ`` times dropped), for its ``--route-top`` leaves with at least
   ``--route-votes``, joined with tier 1's own routed pairs; ``--route tier1`` uses tier 1's
   pairs alone, each leaf expanded to its group.
2. **Chunks along the tree.** Candidate leaves in leaf order (tier 1's merge-tree order, so
   relatives are adjacent and a chunk is a contiguous range) are packed into chunks of at most
   ``--chunk-residues`` residues. Their sequences (each leaf's longest protein) come from the
   index's leaf-ordered store (``--store``, ``store.py``: only the row groups holding them), or
   without one from one streaming pass over ``--members`` that writes them in leaf order. A
   layer group is a single-linkage cluster of the tree, so one contiguous range: the work list
   splits by chunk once, each group in one chunk or two.
3. **Alignment per chunk** (``--method``: ``gapped``, the default, banded
   Smith-Waterman with affine gaps around the anchored seeds, or ``anchored``, ungapped and ~4x
   faster), in batches of at most
   ``--batch-pairs`` pairs, each batch's scores spilled to ``--out``/spill/; then a read-level
   bar at each ``--shares`` where decoy reads called are at most that share of the sample's.
4. **Read classification, ambiguity kept.** A read pair's hits within ``--ambiguous-bits`` of
   its best stay (an ambiguous pair keeps several), shared by EM (``--em-rounds``): weights
   proportional to each leaf's abundance, re-estimated from the weights. Tier 2's presence:
   every leaf with a share in a group with at least ``--min-reads`` expected reads
   (near-identical leaves split a read evenly, which EM cannot break). Routing to an aligner
   downstream comes after, from these calls (``handoff/``).

Memory: one chunk's sequences and seed index (~12 B per residue), one batch of reads
(``--read-batch``; every read pair's frames are written once to ``frames_*.parquet`` and read
back a batch at a time) and one batch of pairs. With ``--route tier1``, the routed reads'
frames come from ``verify.py``'s ``frames_*.parquet`` (recomputed only if absent). The work
list is expanded chunk by chunk. ``--memory-gb`` sets ``--chunk-residues`` and
``--batch-pairs`` from a budget (:func:`budget`), unless they are given.

Writes to ``--out``: ``reads.parquet`` (read, pos, score, weight: each read pair's hits within
``--ambiguous-bits`` of its best and their EM shares), ``cand_verify_<share>.npy`` (tier 2's
presence at each bar), ``verify.json`` (``read_bar_<share>``, as ``verify.py`` names it, so
``recall.py reads`` scores tier 2 as it scores tier 1), ``tier2.json`` (work list, chunks, times,
peak RSS) and, for an
aligner downstream, ``handoff/`` in ``kfp-pileup reference``'s format (:func:`handoff`: the
present leaves as targets, their groups' other leaves as neighbours) with the read pairs
called at the first share's bar (``reads_{1,2}.fq.gz``), for its DIAMOND and ``call``.
"""

import argparse
import gzip
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
import store as leaf_store_mod
import verify
from numpy.typing import NDArray

from kmer_functional_profiler import _core

# BLOSUM62 gapped, open 11 extend 1 (BLAST); verify.LAMBDA, verify.K are the ungapped values
GAPPED_LAMBDA, GAPPED_K = 0.267, 0.041
# bytes per residue of a chunk's sequences, and per pair of a batch (a polars row, the kernel's
# arrays and the scores), measured on the quarter sample with headroom
RESIDUE_BYTES, PAIR_BYTES = 8, 64


def peak_rss_gb() -> float:
    """This process's peak RSS (macOS reports bytes, Linux KiB)."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / (2**30 if sys.platform == "darwin" else 2**20)


def budget(memory_gb: float) -> tuple[int, int]:
    """(chunk residues, batch pairs) for a memory budget: a third each to a chunk's sequences
    and a batch of pairs, the rest left to the routed reads' frames and the work list."""
    b = memory_gb * 2**30 / 3
    return int(b / RESIDUE_BYTES), int(b / PAIR_BYTES)


def groups_of(index: Path) -> NDArray[np.int64]:
    pos_len = np.load(index / "pos_len.npy")
    gfile = index / "layer_group.npy"
    return np.load(gfile).astype(np.int64) if gfile.exists() else np.arange(len(pos_len))


def neighbour_groups(index: Path, group: NDArray[np.int64], d: float) -> pl.DataFrame:
    """(g, nb): pairs of groups joined by a spanning-tree edge at dense-set distance <= d
    (w / (|A| + |B|), |A| ~ length - topology_k + 1), both ways; each group with itself."""
    pos_len = np.load(index / "pos_len.npy").astype(np.int64)
    own = pl.DataFrame({"g": np.unique(group)}).with_columns(nb=pl.col("g"))
    if d <= 0:
        return own
    k_t = int(json.loads((index / "meta.json").read_text())["params"]["topology_k"])
    ea, eb = np.load(index / "edge_a.npy").astype(np.int64), np.load(index / "edge_b.npy")
    eb, w = eb.astype(np.int64), np.load(index / "edge_w.npy").astype(np.float64)
    size = np.maximum(pos_len - k_t + 1, 1)
    ok = (w <= d * (size[ea] + size[eb])) & (group[ea] != group[eb])
    a, b = group[ea[ok]], group[eb[ok]]
    pairs = pl.DataFrame({"g": np.r_[a, b], "nb": np.r_[b, a]})
    return pl.concat([own, pairs]).unique()


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
        .select("protein_id", "pos", "cluster_rep")
    )
    (
        pl.scan_parquet(members)
        .select("protein_id", "sequence")
        .join(prot, on="protein_id")
        .with_columns(n=pl.col("sequence").str.len_bytes())
        .sort("pos", "n", descending=[False, True])
        .unique("pos", keep="first", maintain_order=True)
        .select("pos", "protein_id", "cluster_rep", "sequence")
        .sink_parquet(dst, row_group_size=20_000)
    )


def _masked(s: str) -> bytes:
    return verify.mask(s.encode())


def chunk_seqs(store: Path, lo: int, hi: int) -> dict[int, bytes]:
    """Leaves ``lo..hi``'s sequences from the leaf-ordered store, masked."""
    got = pl.scan_parquet(store).filter(pl.col("pos").is_between(lo, hi)).collect()
    return dict(
        zip(got["pos"].to_list(), (verify.mask(s.encode()) for s in got["sequence"]), strict=True)
    )


def work_list(
    aligned: pl.DataFrame, group: NDArray[np.int64], floor: float, nb: pl.DataFrame | None = None
) -> pl.DataFrame:
    """(read, g): each read's groups with a tier-1 pair scoring at least ``floor``, and with
    ``nb`` (:func:`neighbour_groups`), their neighbours."""
    kept = aligned.filter(pl.col("score") >= floor)
    wl = kept.select("read", g=pl.Series(group[kept["pos"].to_numpy()])).unique()
    if nb is None:
        return wl
    return wl.join(nb, on="g").select("read", g="nb").unique()


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


def all_frames(reads: list[str], dst: Path, batch: int = 100_000) -> int:
    """Every read pair's twelve frames (both mates' six), masked and joined by stops, to
    ``dst`` (read, frames), streamed in batches; returns the read pairs."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    writer, n, rows = None, 0, []

    def flush() -> None:
        nonlocal writer
        t = pa.table({"read": pa.array([r for r, _ in rows], pa.int64()),
                      "frames": pa.array([f for _, f in rows], pa.binary())})  # fmt: skip
        writer = writer or pq.ParquetWriter(dst, t.schema)
        writer.write_table(t)
        rows.clear()

    with gzip.open(reads[0], "rb") as f, gzip.open(reads[1], "rb") as g:
        for i, (x, y) in enumerate(zip(f, g, strict=True)):
            if i % 4 == 1:
                fr = _core.translate_frames(x.strip()) + _core.translate_frames(y.strip())
                rows.append((n, verify.mask(b"*".join(fr))))
                n += 1
                if len(rows) >= batch:
                    flush()
    if rows or writer is None:
        flush()
    if writer:
        writer.close()
    return n


def dense_pairs(
    frames: pl.DataFrame, seqs: dict[int, bytes], a: argparse.Namespace
) -> pl.DataFrame:
    """(read, pos): each read's top ``--route-top`` leaves of the chunk by seed votes
    (``_core.route_reads``: every protein ``--route-k``-mer of the leaves, ``--route-votes`` on
    one diagonal band), the tier-1 trick of short seeds checked by an alignment after."""
    pos = np.array(sorted(seqs), dtype=np.int64)
    r, leaf, _ = _core.route_reads(
        frames["frames"].to_list(), [seqs[x] for x in pos.tolist()], a.route_k,
        verify.protein_classes(), band=8, min_votes=a.route_votes, top=a.route_top,
        max_occ=a.route_max_occ,
    )  # fmt: skip
    return pl.DataFrame({"read": frames["read"].to_numpy()[r], "pos": pos[leaf]})


def align_chunk(
    pairs: pl.DataFrame,
    frames: dict[int, bytes],
    seqs: dict[int, bytes],
    matrix: NDArray[np.int8],
    batch: int,
    method: str = "anchored",
) -> pl.DataFrame:
    """(read, pos, score) for a chunk's pairs, ``batch`` pairs at a time: each batch passes
    only its own reads' frames and leaves' sequences to the kernel."""
    lam, k = (GAPPED_LAMBDA, GAPPED_K) if method == "gapped" else (verify.LAMBDA, verify.K)
    out = []
    pairs = pairs.sort("pos", "read")  # one sketch per run of a leaf
    for start in range(0, pairs.height, batch):
        b = pairs.slice(start, batch)
        reads, ri = np.unique(b["read"].to_numpy(), return_inverse=True)
        pos, pi = np.unique(b["pos"].to_numpy(), return_inverse=True)
        s = _core.pseudo_align(
            [frames[r] for r in reads], [seqs[p] for p in pos.tolist()],
            ri.astype(np.uint32), pi.astype(np.uint32), method, 4, verify.protein_classes(),
            matrix, sketch=(8, 1, 8), band=8, gap=(8, 1), top=2,
        )  # fmt: skip
        bits = (lam * s.astype(np.float64) - math.log(k)) / math.log(2)
        out.append(b.select("read", "pos").with_columns(score=bits))
    return pl.concat(out) if out else pairs.select("read", "pos").with_columns(score=pl.lit(0.0))


def em_weights(calls: pl.DataFrame, rounds: int) -> pl.DataFrame:
    """``calls`` (read, pos, ...: each read's tied leaves) with ``weight``: each read's share
    to each tied leaf, proportional to the leaf's abundance (expected reads over the sample),
    which the shares re-estimate (EM; from equal shares)."""
    w = calls.with_columns(weight=1.0 / pl.len().over("read"))
    for _ in range(rounds):
        theta = w.group_by("pos").agg(theta=pl.col("weight").sum())
        w = (
            w.drop("weight")
            .join(theta, on="pos")
            .with_columns(weight=pl.col("theta") / pl.col("theta").sum().over("read"))
            .drop("theta")
        )
    return w


NEIGHBOURS = 5  # per group: the reference's neighbours, as kfp-pileup's REFERENCE


def handoff(
    out: Path,
    seqs: pl.DataFrame,
    present: NDArray[np.int64],
    group: NDArray[np.int64],
    reads: list[str],
    calls: pl.DataFrame,
) -> None:
    """``out``/handoff/, in ``kfp-pileup reference``'s format (``ref.fa`` with ids s0, s1, ...;
    ``refs.tsv``; ``molecule.txt``), for its DIAMOND and ``call``: the present leaves'
    proteins as targets, up to :data:`NEIGHBOURS` other leaves of each of their groups as
    neighbours (the group as the component); the called read pairs (``reads_{1,2}.fq.gz``, in
    input order) and their classification (``assignments.parquet``: pair, seq_id, unit, score,
    EM weight; an ambiguous pair has several rows), so routing follows it."""
    dst = out / "handoff"
    dst.mkdir(exist_ok=True)
    seqs = seqs.with_columns(component=pl.Series(group[seqs["pos"].to_numpy()]))
    targets = seqs.filter(pl.col("pos").is_in(pl.Series(present).implode()))
    others = seqs.filter(
        pl.col("component").is_in(targets["component"].implode())
        & ~pl.col("pos").is_in(pl.Series(present).implode())
    ).filter(pl.int_range(pl.len()).over("component") < NEIGHBOURS)
    refs = pl.concat(
        [targets.with_columns(role=pl.lit("target")), others.with_columns(role=pl.lit("neighbour"))]
    ).with_columns(
        unit=pl.col("cluster_rep"),
        rep=pl.col("protein_id") == pl.col("cluster_rep"),
        seq_id=pl.format("s{}", pl.int_range(pl.len())),
        length=pl.col("sequence").str.len_bytes(),
    )
    with open(dst / "ref.fa", "w") as f:
        for sid, seq in refs.select("seq_id", "sequence").iter_rows():
            f.write(f">{sid}\n{seq}\n")
    refs.select(
        "unit", "cluster_rep", "protein_id", "rep", "component", "role", "seq_id", "length", "pos"
    ).write_csv(dst / "refs.tsv", separator="\t")
    (dst / "molecule.txt").write_text("protein\n")
    # the read classification, for routing downstream: each called pair's hits (ambiguous
    # pairs several) on the reference, with their EM shares
    calls.join(refs.select("pos", "seq_id", "unit"), on="pos").select(
        pair="read", seq_id="seq_id", unit="unit", score="score", weight="weight"
    ).sort("pair", "seq_id").write_parquet(dst / "assignments.parquet")
    want = set(calls["read"].to_list())
    for m, path in enumerate(reads, start=1):
        with gzip.open(path, "rt") as f, gzip.open(dst / f"reads_{m}.fq.gz", "wt", 1) as g:
            for i, rec in enumerate(zip(f, f, f, f, strict=False)):
                if i in want:
                    g.write("".join(rec))


def main() -> None:
    p = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    p.add_argument("--index", required=True)
    p.add_argument("--members", help="the index's members (sequences), without --store")
    p.add_argument("--store", help="the index's leaf-ordered store (store.py)")
    p.add_argument("--reads", nargs=2, required=True)
    p.add_argument("--query", required=True, help="query.py and verify.py's --out")
    p.add_argument("--share", default="0.01", help="tier 1's read bar share")
    p.add_argument("--route-below", type=float, default=10.0,
                   help="route tier-1 pairs within this many bits of its read bar")  # fmt: skip
    p.add_argument("--neighbours", type=float, default=0.83,
                   help="also route to spanning-tree neighbour groups this close")  # fmt: skip
    p.add_argument("--method", choices=["anchored", "gapped"], default="gapped")
    p.add_argument("--route", choices=["dense", "tier1"], default="dense",
                   help="dense: every read by seed votes; tier1: its pairs")  # fmt: skip
    p.add_argument("--route-k", type=int, default=5)
    p.add_argument("--route-votes", type=int, default=4)
    p.add_argument("--route-top", type=int, default=16)
    p.add_argument("--route-max-occ", type=int, default=2000)
    p.add_argument("--read-batch", type=int, default=200_000)
    p.add_argument("--ambiguous-bits", type=float, default=3.0,
                   help="a read keeps every hit within this many bits of its best")  # fmt: skip
    p.add_argument("--memory-gb", type=float, default=4.0)
    p.add_argument("--chunk-residues", type=int, help="default: from --memory-gb")
    p.add_argument("--batch-pairs", type=int, help="default: from --memory-gb")
    p.add_argument("--em-rounds", type=int, default=20)
    p.add_argument("--min-reads", type=float, default=0.5, help="expected reads for presence")
    p.add_argument("--shares", default="0.01,0.05,0.1")
    p.add_argument("--no-handoff", action="store_true")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    t0 = time.perf_counter()
    index, qdir, out = Path(a.index), Path(a.query), Path(a.out)
    chunk_residues, batch_pairs = budget(a.memory_gb)
    chunk_residues, batch_pairs = a.chunk_residues or chunk_residues, a.batch_pairs or batch_pairs
    spill = out / "spill"
    spill.mkdir(parents=True, exist_ok=True)
    tier1 = query.Index(a.index)
    matrix = verify.blosum62()
    group = groups_of(index)
    nb = neighbour_groups(index, group, a.neighbours) if a.neighbours > 0 else None
    floor = json.loads((qdir / "verify.json").read_text())[f"read_bar_{a.share}"] - a.route_below
    wls = {
        name: work_list(pl.read_parquet(qdir / f), group, floor, nb)
        for name, f in (("sample", "aligned.parquet"), ("decoy", "aligned_decoy.parquet"))
    }
    # the candidates: tier 1's present groups and the groups its pairs route reads to, with
    # their spanning-tree neighbours
    cand = np.load(qdir / f"cand_verify_{a.share}.npy").astype(np.int64)
    cand_g = pl.DataFrame({"g": np.unique(group[cand])})
    if nb is not None:
        cand_g = cand_g.join(nb, on="g").select(g="nb")
    routed_g = pl.concat([cand_g["g"], *(w["g"] for w in wls.values())]).unique().to_numpy()
    leaves = np.flatnonzero(np.isin(group, routed_g))
    length = np.load(index / "pos_len.npy").astype(np.int64)
    chunk_of = chunks(leaves, length, chunk_residues)
    leaf_chunk = pl.DataFrame({"pos": leaves, "g": group[leaves], "chunk": chunk_of})
    info: dict[str, float | int | str] = {
        "method": a.method,
        "route": a.route,
        "neighbours": a.neighbours,
        "candidate_leaves": len(leaves),
        "chunks": int(chunk_of.max()) + 1 if len(leaves) else 0,
        "chunk_residues": chunk_residues,
        "batch_pairs": batch_pairs,
        "residues": int(length[leaves].sum()),
        "route_floor": floor,
    }
    # sequences: the index's leaf-ordered store (store.py) read by position, else one streaming
    # pass over the members writing the candidates' in leaf order
    if a.store:
        ls = leaf_store_mod.LeafStore(a.store)

        def get(pos: NDArray[np.int64]) -> dict[int, bytes]:
            got = ls.get(pos)
            return dict(zip(got["pos"].to_list(), map(_masked, got["sequence"]), strict=True))

        def ids(pos: NDArray[np.int64]) -> pl.DataFrame:
            return ls.get(pos, ids=True)
    else:
        store = out / "seqs.parquet"
        leaf_store(index, a.members, leaves, store)

        def get(pos: NDArray[np.int64]) -> dict[int, bytes]:
            return chunk_seqs(store, int(pos.min()), int(pos.max()))

        def ids(pos: NDArray[np.int64]) -> pl.DataFrame:
            return pl.read_parquet(store).filter(pl.col("pos").is_in(pl.Series(pos).implode()))

    info["store_s"] = time.perf_counter() - t0
    # chunks holding each group (one, or two where a chunk boundary cuts it)
    g_chunk = leaf_chunk.select("g", "chunk").unique()
    dense = a.route == "dense"
    saved = (qdir / "frames_decoy.parquet").exists() and not dense
    reads = {"sample": a.reads, "decoy": [] if saved else query.reversed_reads(a.reads, out)}
    for name, wl in wls.items():
        t1 = time.perf_counter()
        if dense:  # every read pair's frames, once, read back a batch at a time
            fpath = out / f"frames_{name}.parquet"
            info[f"{name}_read_pairs"] = all_frames(reads[name], fpath)
        else:
            frames = routed_frames(tier1, qdir / f"frames_{name}.parquet", reads[name], wl)
        info[f"{name}_frames_s"] = time.perf_counter() - t1
        t1, n_pairs, route_s = time.perf_counter(), 0, 0.0
        parts = wl.join(g_chunk, on="g").partition_by("chunk", as_dict=True)
        for c in range(info["chunks"]):  # type: ignore[arg-type]
            mine = leaf_chunk.filter(pl.col("chunk") == c)
            seqs = get(mine["pos"].to_numpy())
            part = parts.get((c,), pl.DataFrame(schema={"read": pl.Int64, "g": pl.Int64}))
            t1_pairs = part.join(mine.select("g", "pos"), on="g").select("read", "pos")
            if not dense:
                if t1_pairs.height:
                    got = align_chunk(t1_pairs, frames, seqs, matrix, batch_pairs, a.method)
                    got.write_parquet(spill / f"{name}.{c:05d}.parquet")
                    n_pairs += t1_pairs.height
                continue
            for b, batch in enumerate(pl.read_parquet(fpath).iter_slices(a.read_batch)):
                t2 = time.perf_counter()
                pairs = dense_pairs(batch, seqs, a)
                route_s += time.perf_counter() - t2
                lo, hi = int(batch["read"][0]), int(batch["read"][-1])
                mine_t1 = t1_pairs.filter(pl.col("read").is_between(lo, hi))
                pairs = pl.concat([pairs, mine_t1.cast(pairs.schema)]).unique()
                if not pairs.height:
                    continue
                fr = dict(zip(batch["read"].to_list(), batch["frames"].to_list(), strict=True))
                got = align_chunk(pairs, fr, seqs, matrix, batch_pairs, a.method)
                got.write_parquet(spill / f"{name}.{c:05d}.{b:05d}.parquet")
                n_pairs += pairs.height
        info |= {f"{name}_pairs": n_pairs, f"{name}_route_s": route_s,
                 f"{name}_align_s": time.perf_counter() - t1 - route_s}  # fmt: skip
        if dense:
            fpath.unlink()
        else:
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
        # a read's hits within --ambiguous-bits of its best: ambiguous reads keep several
        best[name] = lf.filter(
            pl.col("score") >= pl.col("score").max().over("read") - a.ambiguous_bits
        ).collect()
    per_read = {n: b.group_by("read").agg(pl.col("score").max()) for n, b in best.items()}
    bars: dict[str, float] = {}
    for s in a.shares.split(","):
        bars[f"read_bar_{float(s)}"] = query.decoy_bar(per_read["sample"], per_read["decoy"],
                                                       float(s))[0]  # fmt: skip
    # EM over the reads called at the loosest bar; presence per bar from the shares
    loosest = min(bars.values())
    top = pl.col("score").max().over("read")
    calls = em_weights(best["sample"].filter(top >= loosest), a.em_rounds)
    calls.sort("read", "pos").write_parquet(out / "reads.parquet")
    first = None
    for s in a.shares.split(","):
        # a group is present with at least min_reads expected reads (near-identical leaves split
        # a read's weight evenly, which EM cannot break), and so is each of its leaves with any
        at = calls.filter(pl.col("score").max().over("read") >= bars[f"read_bar_{float(s)}"])
        expected = at.group_by("pos").agg(pl.col("weight").sum())
        expected = expected.with_columns(g=pl.Series(group[expected["pos"].to_numpy()]))
        ok = expected.group_by("g").agg(n=pl.col("weight").sum()).filter(pl.col("n") >= a.min_reads)
        present = np.sort(expected.filter(pl.col("g").is_in(ok["g"].implode()))["pos"].to_numpy())
        np.save(out / f"cand_verify_{float(s)}.npy", present)
        if first is None:
            first = (present, at)
    (out / "verify.json").write_text(json.dumps(bars, indent=1) + "\n")
    if not a.no_handoff and first is not None:
        t1 = time.perf_counter()
        handoff(out, ids(leaves), first[0], group, a.reads, first[1])
        info["handoff_s"] = time.perf_counter() - t1
    shutil.rmtree(spill)
    if not a.store:
        store.unlink()
    for f in out.glob("decoy_*.fq.gz"):
        f.unlink()
    info |= bars | {"total_s": time.perf_counter() - t0, "peak_rss_gb": peak_rss_gb()}
    (out / "tier2.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info))


if __name__ == "__main__":
    main()
