"""Tree index, phase 2: counting (docs/tree-index-plan.md, Querying a metagenome, steps 1-2).

- **Lookup:** the sample's k-mers (reads in translated frames as today's query, or protein
  sequences) with hash <= t\\_cap, or equal to one of the few keys beyond it (top-ups of short
  proteins), counted per distinct hash and probed against the sorted key array.
- **Counting:** each hit k-mer's postings add their sign over their range. A sweep over the
  sorted range endpoints gives piecewise-constant counts over leaf order; internal nodes take
  the postings whose range contains theirs (``_core.contained_sums``). Two counts each: distinct
  hit k-mers (``hit_kmers``) and their occurrences in the sample (``hits``).

Writes to ``--out``: ``leaves.parquet`` (pos, size, hit\\_kmers, hits, containment) and
``nodes.parquet`` (node, lo, hi, level, size, hit\\_kmers, hits, containment), rows with a hit
k-mer only; ``hits.parquet`` (hash, n: the sample's k-mers that are keys, for ``report.py``); and
``query.json`` (counts and timings).
"""

import argparse
import json
import time
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from numpy.typing import NDArray

from kmer_functional_profiler import _core


class Index:
    """A ``tree-v0`` index, memory-mapped."""

    def __init__(self, path: str | Path) -> None:
        path = Path(path)
        self.meta = json.loads((path / "meta.json").read_text())
        load = lambda name: np.load(path / f"{name}.npy", mmap_mode="r")  # noqa: E731
        self.hash, self.lo, self.hi, self.sign = (load(x) for x in ("hash", "lo", "hi", "sign"))
        self.pos_size, self.node_size = load("pos_size"), load("node_size")
        self.node_lo, self.node_hi, self.node_level = (
            load(x) for x in ("node_lo", "node_hi", "node_level")
        )
        # keys are sampled up to t_cap, except a few top-ups beyond it (short proteins)
        self.sketch_max = _core.max_hash(float(self.meta["params"]["t_cap"]))
        tail = np.asarray(self.hash[np.searchsorted(self.hash, self.sketch_max, "right") :])
        self.beyond = np.unique(tail)

    @property
    def k(self) -> int:
        return int(self.meta["params"]["k"])

    @property
    def max_hash(self) -> int:
        return int(self.hash[-1]) if len(self.hash) else 0

    def sampled(self, h: NDArray[np.uint64]) -> NDArray[np.uint64]:
        """The hashes a query keeps: up to t\\_cap, or a key beyond it."""
        at = np.minimum(np.searchsorted(self.beyond, h), max(len(self.beyond) - 1, 0))
        hit = self.beyond[at] == h if len(self.beyond) else np.zeros(len(h), dtype=bool)
        return h[(h <= self.sketch_max) | hit]


def sample_counts(batches: Iterable[NDArray[np.uint64]]) -> pl.DataFrame:
    """(hash, n): occurrences of each distinct hash over the batches."""
    parts = [pl.DataFrame({"hash": b}).group_by("hash").len("n") for b in batches]
    if not parts:
        return pl.DataFrame(schema={"hash": pl.UInt64, "n": pl.UInt32})
    return pl.concat(parts).group_by("hash").agg(pl.col("n").sum()).sort("hash")


def read_hashes(index: Index, r1: str, r2: str | None, frames: str) -> Iterable[NDArray[Any]]:
    for b in _core.FastxHits(r1, r2, k=index.k, frames=frames, max_hash=index.max_hash):
        yield index.sampled(b["hash"])


def protein_hashes(index: Index, fasta: str) -> Iterable[NDArray[Any]]:
    seqs, cur = [], []
    for line in Path(fasta).read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur).encode())
            cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur).encode())
    # one sequence at a time: hash_proteins dedups within a call, so n counts sequences
    for s in seqs:
        yield index.sampled(_core.hash_proteins([s], index.k, max_hash=index.max_hash)["hash"])


def count(index: Index, sample: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame, dict]:
    """Hit counts per leaf position and per internal node, for a (hash, n) sample."""
    h = sample["hash"].to_numpy()
    first = np.searchsorted(index.hash, h, "left")
    last = np.searchsorted(index.hash, h, "right")
    per = last - first
    rows = np.repeat(first - np.cumsum(per) + per, per) + np.arange(per.sum())
    n = np.repeat(sample["n"].to_numpy().astype(np.int64), per)
    sign = np.asarray(index.sign[rows], dtype=np.int64)
    lo, hi = np.asarray(index.lo[rows]), np.asarray(index.hi[rows])
    w = {"hit_kmers": sign, "hits": sign * n}
    # leaves: one sweep over the endpoints, then the positive segments expanded
    ends, at = np.unique(np.concatenate([lo, hi]), return_inverse=True)
    seg = {}
    for name, x in w.items():
        delta = np.zeros(len(ends), dtype=np.int64)
        np.add.at(delta, at, np.concatenate([x, -x]))
        seg[name] = np.cumsum(delta)[:-1] if len(ends) else delta
    on = np.flatnonzero(seg["hit_kmers"] > 0)
    lens = (ends[on + 1] - ends[on]).astype(np.int64)
    pos = np.repeat(ends[on].astype(np.int64) - np.cumsum(lens) + lens, lens) + np.arange(
        lens.sum()
    )
    leaves = pl.DataFrame(
        {
            "pos": pos,
            "size": np.asarray(index.pos_size[pos], dtype=np.int64),
            "hit_kmers": np.repeat(seg["hit_kmers"][on], lens),
            "hits": np.repeat(seg["hits"][on], lens),
        }
    )
    # ponytail: every internal node is a query, O(N log P); touched ancestors only at scale
    lo32, hi32 = lo.astype(np.uint32), hi.astype(np.uint32)
    node_hit = {
        name: _core.contained_sums(lo32, hi32, x, index.node_lo, index.node_hi)
        for name, x in w.items()
    }
    keep = np.flatnonzero(node_hit["hit_kmers"] > 0)
    nodes = pl.DataFrame(
        {
            "node": keep,
            "lo": np.asarray(index.node_lo[keep]),
            "hi": np.asarray(index.node_hi[keep]),
            "level": np.asarray(index.node_level[keep]),
            "size": np.asarray(index.node_size[keep], dtype=np.int64),
            "hit_kmers": node_hit["hit_kmers"][keep],
            "hits": node_hit["hits"][keep],
        }
    )
    leaves, nodes = (
        t.with_columns(containment=pl.col("hit_kmers") / pl.col("size")) for t in (leaves, nodes)
    )
    stats = {
        "sampled_hashes": len(h),
        "hit_hashes": int((per > 0).sum()),
        "hit_postings": len(rows),
        "leaves_hit": leaves.height,
        "nodes_hit": nodes.height,
    }
    return leaves, nodes, stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--index", required=True)
    parser.add_argument("--out", required=True)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--reads", nargs="+", help="R1 [R2] FASTQ")
    src.add_argument("--proteins", help="protein FASTA")
    parser.add_argument("--frames", default="stopfree")
    a = parser.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    index = Index(a.index)
    t0 = time.perf_counter()
    if a.reads:
        batches = read_hashes(index, a.reads[0], (a.reads[1:] or [None])[0], a.frames)
    else:
        batches = protein_hashes(index, a.proteins)
    sample = sample_counts(batches)
    lookup_s = time.perf_counter() - t0
    leaves, nodes, stats = count(index, sample)
    leaves.write_parquet(out / "leaves.parquet")
    nodes.write_parquet(out / "nodes.parquet")
    h = sample["hash"].to_numpy()
    at = np.minimum(np.searchsorted(index.hash, h), max(len(index.hash) - 1, 0))
    keyed = index.hash[at] == h if len(index.hash) else np.zeros(len(h), dtype=bool)
    sample.filter(pl.Series(keyed)).write_parquet(out / "hits.parquet")
    stats |= {"lookup_s": lookup_s, "total_s": time.perf_counter() - t0}
    (out / "query.json").write_text(json.dumps(stats, indent=1) + "\n")
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
