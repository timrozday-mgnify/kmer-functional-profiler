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
``query.json`` (counts and timings). With reads and an index with key bits:
``seeds.parquet`` (pos, score: the best evidence any one read pair gives each leaf, in bits;
:func:`evidence`), ``evidence.parquet`` (the curve :func:`choose_bits` reads) and the chosen
``evidence_bits`` in ``query.json``.
"""

import argparse
import gzip
import json
import time
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Final

import alphabets
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
        if (path / "keys.npy").exists():
            self.keys, self.key_bits, self.key_k = (
                np.load(path / f"{x}.npy") for x in ("keys", "key_bits", "key_k")
            )

    @property
    def k(self) -> int:
        return int(self.meta["params"]["k"])

    def keyed(self, h: NDArray[np.uint64]) -> NDArray[np.uint64]:
        """The hashes that are keys."""
        at = np.minimum(np.searchsorted(self.hash, h), max(len(self.hash) - 1, 0))
        return h[np.asarray(self.hash[at]) == h] if len(self.hash) else h[:0]

    @property
    def max_hash(self) -> int:
        return int(self.hash[-1]) if len(self.hash) else 0

    def sampled(self, h: NDArray[np.uint64]) -> NDArray[np.uint64]:
        """The hashes a query keeps: up to t\\_cap, or a key beyond it."""
        at = np.minimum(np.searchsorted(self.beyond, h), max(len(self.beyond) - 1, 0))
        hit = self.beyond[at] == h if len(self.beyond) else np.zeros(len(h), dtype=bool)
        return h[(h <= self.sketch_max) | hit]


def read_pairs(
    index: Index,
    r1: str,
    r2: str | None,
    frames: str,
    ladder_frames: str | None = None,
    mode: str = "recall",
) -> tuple[pl.DataFrame, dict[int, int], int]:
    """Distinct (read, hash) of the sample's k-mers that are keys, at the index's k and each
    ladder k; the read windows hashed per k; and the read pairs."""
    parts, windows, n_reads = [], {}, 0
    todo = [(index.k, "protein", 0)] + [
        (k, alpha, code) for k, _, alpha, code in steps(index, mode)
    ]
    for k, alpha, code in todo:
        f = (ladder_frames or frames) if code else frames
        it = _core.FastxHits(r1, r2, k=k, frames=f, alphabet=alpha)
        for b in it:
            h = salted(b["hash"], code) if code else b["hash"]
            key = code or k
            windows[key] = windows.get(key, 0) + len(h)
            hit = np.isin(h, index.keyed(h))
            parts.append(pl.DataFrame({"read": b["read"][hit], "hash": h[hit]}))
        n_reads = max(n_reads, it.n_reads)
    return pl.concat(parts).unique(), windows, n_reads


def evidence_chunked(
    index: Index, pairs: pl.DataFrame, log2_windows: float, clamp: bool, chunk: int = 20_000
) -> pl.DataFrame:
    """:func:`evidence` over read pairs ``chunk`` at a time (memory bounded; the same result,
    since evidence is per read pair)."""
    parts = []
    reads = pairs["read"].to_numpy()
    for lo in range(0, int(reads.max(initial=0)) + 1, chunk):
        sub = pairs.filter((pl.col("read") >= lo) & (pl.col("read") < lo + chunk))
        if sub.height:
            parts.append(evidence(index, sub, log2_windows, floor=0.0, clamp=clamp))
    if not parts:
        return pl.DataFrame(schema={"pos": pl.Int64, "score": pl.Float64})
    return pl.concat(parts).group_by("pos").agg(pl.col("score").max()).sort("pos")


def evidence(
    index: Index, pairs: pl.DataFrame, log2_windows: float, floor: float, clamp: bool = False
) -> pl.DataFrame:
    """(pos, score): each leaf position's best evidence from one read pair, the read pair's
    hit k-mers on it summed in bits less log2(windows per read pair) per hit beyond the
    first (each further hit had that many places to fall by chance); scores >= ``floor``.
    ``clamp``: a hit's bits are at least 0 (common k-mers add nothing rather than subtract)."""
    h = pairs["hash"].to_numpy()
    bits = index.key_bits[np.searchsorted(index.keys, h)].astype(np.float64)
    if clamp:
        bits = np.maximum(bits, 0.0)
    first = np.searchsorted(index.hash, h, "left")
    per = np.searchsorted(index.hash, h, "right") - first
    rows = np.repeat(first - np.cumsum(per) + per, per) + np.arange(per.sum())
    read = np.repeat(pairs["read"].to_numpy(), per)
    sign = np.asarray(index.sign[rows], dtype=np.float64)
    b = np.repeat(bits, per) * sign
    x = np.concatenate([np.asarray(index.lo[rows]), np.asarray(index.hi[rows])]).astype(np.int64)
    r = np.concatenate([read, read])
    o = np.lexsort((x, r))
    x, r = x[o], r[o]
    cb = np.cumsum(np.concatenate([b, -b])[o])  # signs balance per k-mer: sums restart per read
    cm = np.cumsum(np.concatenate([sign, -sign])[o])
    score = cb - np.maximum(cm - 1, 0) * log2_windows
    on = np.flatnonzero(
        (cm[:-1] > 0.5) & (r[1:] == r[:-1]) & (x[1:] > x[:-1]) & (score[:-1] >= floor)
    )
    n = x[on + 1] - x[on]
    pos = np.repeat(x[on] - np.cumsum(n) + n, n) + np.arange(n.sum())
    seg = pl.DataFrame({"pos": pos, "score": np.repeat(score[on], n)})
    return seg.group_by("pos").agg(pl.col("score").max()).sort("pos")


def reversed_reads(paths: list[str], out: Path) -> list[str]:
    """Decoy reads: each read reversed (not complemented), keeping its composition and low
    complexity but not its coding homology."""
    made = []
    for i, path in enumerate(paths):
        dst = out / f"decoy_{i + 1}.fq.gz"
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as f, gzip.open(dst, "wt", compresslevel=1) as g:
            for name, seq, plus, qual in zip(f, f, f, f, strict=False):
                g.write(f"{name}{seq.rstrip()[::-1]}\n{plus}{qual.rstrip()[::-1]}\n")
        made.append(str(dst))
    return made


def decoy_bar(
    seeds: pl.DataFrame, decoys: pl.DataFrame, share: float
) -> tuple[float, list[tuple[float, float, int]]]:
    """The least score at which leaves seeded by decoy reads are at most ``share`` of the
    leaves seeded by the sample (target-decoy); and the curve (score, decoy seeds, seeds)."""
    sc, dc = np.sort(seeds["score"].to_numpy()), np.sort(decoys["score"].to_numpy())
    curve = []
    for bar in np.arange(0.0, 160.0, 0.25):
        n, d = len(sc) - np.searchsorted(sc, bar), len(dc) - np.searchsorted(dc, bar)
        curve.append((float(bar), float(d), int(n)))
    ok = [c for c in curve if c[1] <= share * max(c[2], 1)]
    return (ok[0][0] if ok else curve[-1][0]), curve


def choose_bits(
    index: Index, seeds: pl.DataFrame, windows: dict[int, int], share: float
) -> tuple[float, list[tuple[float, float, int]]]:
    """The least evidence B at which the expected chance seeds, sum over keys with bits >= B
    of (read windows of its k) x 2^-bits, are at most ``share`` of the seeds scoring >= B;
    and the curve (B, expected chance seeds, seeds). ponytail: chance single hits only;
    multi-hit chance seeds are fewer by the per-hit window charge."""
    w = np.array([windows.get(int(k), 0) for k in index.key_k], dtype=np.float64)
    w[index.key_k == 0] = windows.get(index.k, 0)
    o = np.argsort(-index.key_bits.astype(np.float64))
    kb = index.key_bits[o].astype(np.float64)
    cum = np.cumsum((w[o] * 2.0**-kb)[::-1])[::-1]  # chance seeds from keys with bits >= kb
    sc = np.sort(seeds["score"].to_numpy())
    curve = []
    for bar in np.arange(np.floor(sc.min()) if len(sc) else 0, 120.0, 0.25):
        at = np.searchsorted(-kb, -bar, "right")  # keys with bits >= bar: the first `at`
        chance = float(cum[0] - (cum[at] if at < len(cum) else 0.0)) if at else 0.0
        n = len(sc) - np.searchsorted(sc, bar)
        curve.append((float(bar), chance, int(n)))
    ok = [c for c in curve if c[1] <= share * max(c[2], 1)]
    return (ok[0][0] if ok else curve[-1][0] if curve else 0.0), curve


ALPHABETS: Final = ("protein", "murphy10", "dayhoff", *alphabets.REDUCED)


def ladder(spec: str) -> list[tuple[int, int, str, int]]:
    """``build.py --ladder``'s steps, ``[alphabet:]k:window``: (k, window, alphabet, code),
    code = k + 32 x the alphabet's place in ALPHABETS, which keys the step's hashes and
    windows (a protein step's code is its k)."""
    out = []
    for step in filter(None, spec.split(",")):
        *alpha, k, w = [x for x in step.split(":") if x not in FLAGS]
        name = alpha[0] if alpha else "protein"
        out.append((int(k), int(w), name, int(k) + 32 * ALPHABETS.index(name)))
    return out


FLAGS: Final = ("rare",)


def ladder_flags(spec: str) -> list[set[str]]:
    """Per step, its flags: ``rare`` (each window's rarest k-mer in the database, not its
    least hash)."""
    return [set(step.split(":")) & set(FLAGS) for step in filter(None, spec.split(","))]


def steps(index: Index, mode: str = "recall") -> list[tuple[int, int, str, int]]:
    """The index's ladder and layer steps a query uses: all of them in ``recall`` mode; in
    ``precision`` mode the ladder only (layers are not yet tested for precision or F1)."""
    params = index.meta["params"]
    layers = ladder(params.get("layers", "") or "") if mode == "recall" else []
    return ladder(params.get("ladder", "")) + layers


def salted(h: NDArray[np.uint64], k: int) -> NDArray[np.uint64]:
    """Hashes of a ladder step moved to their own domain (a bijection per step code ``k``):
    k-mer hashes do not depend on k or alphabet, so a 9-mer can share an 11-mer's hash."""
    with np.errstate(over="ignore"):
        x = (h ^ np.uint64(k * 0x9E3779B97F4A7C15 % 2**64)) * np.uint64(0xBF58476D1CE4E5B9)
    return x ^ (x >> np.uint64(31))


def sample_counts(batches: Iterable[NDArray[np.uint64]]) -> pl.DataFrame:
    """(hash, n): occurrences of each distinct hash over the batches."""
    parts = [pl.DataFrame({"hash": b}).group_by("hash").len("n") for b in batches]
    if not parts:
        return pl.DataFrame(schema={"hash": pl.UInt64, "n": pl.UInt32})
    return pl.concat(parts).group_by("hash").agg(pl.col("n").sum()).sort("hash")


def read_hashes(
    index: Index, r1: str, r2: str | None, frames: str, mode: str = "recall"
) -> Iterable[NDArray[Any]]:
    for b in _core.FastxHits(r1, r2, k=index.k, frames=frames, max_hash=index.max_hash):
        yield index.sampled(b["hash"])
    for k, _, alpha, code in steps(index, mode):  # every k-mer that is a key
        for b in _core.FastxHits(r1, r2, k=k, frames=frames, alphabet=alpha):
            yield index.keyed(salted(b["hash"], code))


def protein_hashes(index: Index, fasta: str, mode: str = "recall") -> Iterable[NDArray[Any]]:
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
        h = index.sampled(_core.hash_proteins([s], index.k, max_hash=index.max_hash)["hash"])
        for k, _, alpha, code in steps(index, mode):
            g = _core.hash_proteins([s], k, alphabet=alpha)["hash"]
            h = np.concatenate([h, index.keyed(salted(g, code))])
        yield h


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
    parser.add_argument("--frames", default="all")
    parser.add_argument(
        "--chance-share", type=float, default=0.01,
        help="reads: seeds by evidence, at the least bits whose expected chance seeds are at "
        "most this share of the seeds (seeds.parquet)",
    )  # fmt: skip
    parser.add_argument("--ladder-frames", help="frames for ladder k-mers (default: --frames)")
    parser.add_argument(
        "--mode", choices=["recall", "precision"], default="recall",
        help="recall: every layer of the index; precision: the ladder only, no layers "
        "(layers are not yet tested for precision or F1)",
    )  # fmt: skip
    parser.add_argument(
        "--decoy", action="store_true",
        help="count chance seeds on the reads reversed (target-decoy) instead of modelling them",
    )  # fmt: skip
    a = parser.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    index = Index(a.index)
    t0 = time.perf_counter()
    if a.reads:
        batches = read_hashes(index, a.reads[0], (a.reads[1:] or [None])[0], a.frames, a.mode)
    else:
        batches = protein_hashes(index, a.proteins, a.mode)
    sample = sample_counts(batches)
    lookup_s = time.perf_counter() - t0
    leaves, nodes, stats = count(index, sample)
    leaves.write_parquet(out / "leaves.parquet")
    nodes.write_parquet(out / "nodes.parquet")
    h = sample["hash"].to_numpy()
    at = np.minimum(np.searchsorted(index.hash, h), max(len(index.hash) - 1, 0))
    keyed = index.hash[at] == h if len(index.hash) else np.zeros(len(h), dtype=bool)
    sample.filter(pl.Series(keyed)).write_parquet(out / "hits.parquet")
    if a.reads and hasattr(index, "keys"):
        t1 = time.perf_counter()
        pairs, windows, n_reads = read_pairs(
            index, a.reads[0], (a.reads[1:] or [None])[0], a.frames, a.ladder_frames, a.mode
        )
        log2_w = float(np.log2(max(windows.get(index.k, 1) / max(n_reads, 1), 2)))
        seeds = evidence_chunked(index, pairs, log2_w, clamp=a.decoy)
        if a.decoy:  # the bar from reversed reads: chance seeds counted, not modelled
            dr = reversed_reads(a.reads, out)
            d_pairs, _, _ = read_pairs(index, *dr, a.frames, a.ladder_frames, a.mode)
            decoys = evidence_chunked(index, d_pairs, log2_w, clamp=True)
            decoys.write_parquet(out / "decoy_seeds.parquet")
            bar, curve = decoy_bar(seeds, decoys, a.chance_share)
            stats["decoy_seeds"] = int((decoys["score"] >= bar).sum())
        else:
            bar, curve = choose_bits(index, seeds, windows, a.chance_share)
        seeds.write_parquet(out / "seeds.parquet")
        pl.DataFrame(curve, schema=["bits", "chance", "seeds"], orient="row").write_parquet(
            out / "evidence.parquet"
        )
        stats |= {"evidence_bits": bar, "seeds": int((seeds["score"] >= bar).sum()),
                  "read_pairs": n_reads, "windows": windows,
                  "evidence_s": time.perf_counter() - t1}  # fmt: skip
    stats |= {"lookup_s": lookup_s, "total_s": time.perf_counter() - t0}
    (out / "query.json").write_text(json.dumps(stats, indent=1) + "\n")
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
