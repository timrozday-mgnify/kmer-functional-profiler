"""Tree index, phase 4b: a per-sample table along tier 1's tree against the flat table it would
replace (docs/tree-index-plan.md, Uses beyond tier 1), on tier 2's reference (``handoff/``).

Per ``--k``: every k-mer of the reference's leaves, as today's flat table (``kfp_pileup.route``'s
``kmer_table``, its time) and as Fitch's changes on the merge tree restricted to those leaves
(no spanning tree, nothing built but the Fitch pass). Then the changes are swept back to
(k-mer, leaf), which must equal the flat rows exactly: routing's per-block table and a
per-leaf table for reclassification both follow from it. The reads' k-mers (all frames) give the
postings each table would visit. Writes ``--out`` (TSV, one row per k).
"""

import argparse
import tempfile
import time
from pathlib import Path

import numpy as np
import polars as pl
from numpy.typing import NDArray

from kfp_pileup.route import blocks, kmer_table
from kmer_functional_profiler import _core
from kmer_functional_profiler.index import read_fasta

FLAT_BYTES = 8  # a packed table's entry (route.py's kmer_table)
TREE_BYTES = 12  # ponytail: a packed key's 4 B and lo, hi as u32, the sign in a bit; not measured


def restricted_forest(
    pos: NDArray[np.int64], node_lo: NDArray, node_hi: NDArray
) -> tuple[NDArray[np.uint32], NDArray[np.uint32], NDArray[np.uint32], NDArray[np.int64]]:
    """The merge forest over the sorted leaf positions ``pos``: (left, right, tops, and the
    leaf range per node, ``ranges[v]`` = lo, hi). A node keeps its range over ``pos``; one with
    a single child left, whose range repeats that child's, goes."""
    lo, hi = np.searchsorted(pos, node_lo), np.searchsorted(pos, node_hi)
    keep = hi - lo >= 2
    got = np.unique(np.c_[lo[keep], hi[keep]], axis=0)
    got = got[np.argsort(got[:, 1] - got[:, 0], kind="stable")]  # children before parents
    n = len(pos)
    top, end = list(range(n)), [*range(1, n + 1), *got[:, 1].tolist()]
    left, right = np.zeros(len(got), np.uint32), np.zeros(len(got), np.uint32)

    def walk(a: int, b: int) -> list[int]:
        out = []
        while a < b:
            out.append(top[a])
            a = end[top[a]]
        return out

    for i, (a, b) in enumerate(got.tolist()):
        kids = walk(a, b)
        if len(kids) != 2:
            raise ValueError(f"node ({a}, {b}) has {len(kids)} children: not a laminar binary tree")
        left[i], right[i] = kids
        top[a] = n + i
    ranges = np.r_[np.c_[np.arange(n), np.arange(1, n + 1)], got].astype(np.int64)
    return left, right, np.array(walk(0, n), dtype=np.uint32), ranges


def sweep(
    h: NDArray[np.uint64], lo: NDArray[np.int64], hi: NDArray[np.int64], sign: NDArray[np.int8]
) -> tuple[NDArray[np.uint64], NDArray[np.int64]]:
    """(k-mer, leaf) rows from a k-mer's nested gains and losses: a leaf holds it where the
    changes covering the leaf sum above zero."""
    eh = np.r_[h, h]
    at = np.r_[lo, hi]
    d = np.r_[sign, -sign].astype(np.int64)
    o = np.lexsort((at, eh))
    eh, at, depth = eh[o], at[o], np.cumsum(d[o])
    same = eh[1:] == eh[:-1]
    start = at[:-1][same & (depth[:-1] > 0)]
    stop = at[1:][same & (depth[:-1] > 0)]
    which = eh[:-1][same & (depth[:-1] > 0)]
    n = stop - start
    leaf = np.repeat(start, n) + np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)
    return np.repeat(which, n), leaf


def measure(a: argparse.Namespace) -> None:
    ref, index = Path(a.ref), Path(a.index)
    refs = pl.read_csv(ref / "refs.tsv", separator="\t").sort("pos")
    found = dict(read_fasta(ref / "ref.fa"))
    seqs = [found[s].encode() for s in refs["seq_id"]]
    pos = refs["pos"].to_numpy()
    block = blocks(refs, a.block_letters)
    rows = []
    for k in a.k:
        t = time.perf_counter()
        with tempfile.TemporaryDirectory() as scratch:
            kmer_table(seqs, block, k, "protein", Path(scratch))
        flat_s = time.perf_counter() - t
        t = time.perf_counter()
        got = _core.hash_proteins(seqs, k)
        kmers = pl.DataFrame({"leaf": got["seq"], "h": got["hash"]}).unique().sort("leaf", "h")
        hash_s = time.perf_counter() - t
        leaf, h = kmers["leaf"].to_numpy().astype(np.int64), kmers["h"].to_numpy()
        n = len(seqs)
        t = time.perf_counter()
        left, right, tops, ranges = restricted_forest(
            pos, np.load(index / "node_lo.npy"), np.load(index / "node_hi.npy")
        )
        offsets = np.r_[0, np.cumsum(np.bincount(leaf, minlength=n))].astype(np.uint64)
        ch = _core.tree_schemes(
            n, left, right, tops, offsets, h, np.zeros(n, np.uint8),
            np.array([2**64 - 1], np.uint64), emit=True,
        )  # fmt: skip
        tree_s = time.perf_counter() - t
        t = time.perf_counter()
        node = ch["change_node"].astype(np.int64)
        sign = np.where(ch["change_gain"], 1, -1).astype(np.int8)
        back_h, back_leaf = sweep(ch["change_hash"], ranges[node, 0], ranges[node, 1], sign)
        sweep_s = time.perf_counter() - t
        back = pl.DataFrame({"leaf": back_leaf, "h": back_h}).sort("leaf", "h")
        if not back.equals(kmers.with_columns(pl.col("leaf").cast(pl.Int64))):
            raise ValueError(f"k {k}: the sweep's rows differ from the flat table's")
        by_block = kmers.with_columns(block=pl.Series(block)[kmers["leaf"]]).unique(["h", "block"])
        # postings visited by the reads' k-mers, per k-mer occurrence
        n_flat = kmers.group_by("h").agg(f=pl.len())
        n_tree = pl.DataFrame({"h": ch["change_hash"]}).group_by("h").agg(c=pl.len())
        n_block = by_block.group_by("h").agg(b=pl.len())
        hits = _core.FastxHits(*a.reads, k=k, frames="all", alphabet="protein")
        reads = (
            pl.DataFrame({"h": np.concatenate([b["hash"] for b in hits])})
            .join(n_flat, on="h")
            .join(n_tree, on="h")
            .join(n_block, on="h")
        )
        rows.append(
            {
                "k": k,
                "leaves": n,
                "residues": sum(map(len, seqs)),
                "flat_entries": kmers.height,
                "block_entries": by_block.height,
                "tree_changes": len(node),
                "x_entries": kmers.height / max(len(node), 1),
                "flat_gb": FLAT_BYTES * kmers.height / 1e9,
                "tree_gb": TREE_BYTES * len(node) / 1e9,
                "route_table_s": flat_s,
                "hash_s": hash_s,
                "fitch_s": tree_s,
                "sweep_s": sweep_s,
                "read_hits": reads.height,
                "visit_flat": int(reads["f"].sum()),
                "visit_block": int(reads["b"].sum()),
                "visit_tree": int(reads["c"].sum()),
            }  # fmt: skip
        )
        print(rows[-1], flush=True)
    pl.DataFrame(rows).write_csv(a.out, separator="\t")


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--ref", required=True, help="tier2.py's handoff/")
    parser.add_argument("--index", required=True)
    parser.add_argument("--reads", nargs=2, required=True)
    parser.add_argument("--k", type=int, nargs="+", default=[8])
    parser.add_argument("--block-letters", type=int, default=2_000_000, help="as kfp-pileup route")
    parser.add_argument("--out", required=True)
    measure(parser.parse_args())


if __name__ == "__main__":
    main()
