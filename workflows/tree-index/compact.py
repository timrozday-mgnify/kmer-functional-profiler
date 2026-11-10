"""Compact route index: MGnify90 units, k-mer keys on each unit's representative, nothing kept
per protein (docs/tree-index-plan.md, Progress log: the compact route, 2026-11-10).

**Sequences** (``--seqs``, built once and shared by every index size, ``seqs``):

- ``<shard>.fa``: the shard's proteins as FASTA, unit by unit, the unit's representative
  (MGnify's ``cluster_rep``, else the longest member) first. MGnify's proteins reordered: the
  sequences, not the index. Units go to shards by a hash of the representative.
- ``<shard>.units.npy`` (u64): each unit's first byte in that file, and the end;
  ``<shard>.members.npy`` (u32): its members.

**Index** (``--index``, one per size: one key per ``--rate`` representative residues, 32 the
standard, 16 and 8 the larger):

- ``shards.npy`` (i64): each shard's first global unit, and the units in all.
- ``keys/``: one posting per (key, unit), sorted by the key's code (its k-mer's base-20 rank,
  :func:`codes`): ``bucket.npy`` (u64: where each bucket of ``code >> suffix_bits`` starts),
  ``suffix.npy`` (the code's low bits: u16 at release), ``unit.npy`` (u32) and ``bits.npy``
  (u8 quarter-bits: -log2 q0 - log2 g, g the units holding the key). 7 B per posting at
  release (2^27 buckets, 12-bit suffixes); ``meta.json`` (``seqs``: the sequences' directory,
  relative to the index).

**Build:** ``seqs`` per shard, once; then per size, ``shard`` per shard (independent): the keys,
each representative's k-mer window minimizers (window ``--window``, each window's rarest k-mer
among the shard's representatives) thinned to sqrt(members) x a mean set so the
representatives carry one key per ``--rate`` residues (``build.allocate``); and ``merge``: the
shards' keys by code partition, each key's units counted for its bits.

**Query:** tier 1 (each read pair's 12 frames' k-mers looked up; a (read pair, unit) with at
least ``--floor`` bits: its keys' bits summed less log2 windows per hit beyond the first, as
``query.evidence``), then confirmation (``_core.pseudo_align`` of those pairs against the unit's
representative; a unit is confirmed when its best pair beats the decoys' bar, the reads
reversed, at ``--share`` and ``--chance-scale``). Writes ``units.tsv`` (unit, score: its best
pair, protein_id: its representative), and for recruitment, a step of its own (``recruit.py``):
``confirmed.fa`` (the confirmed representatives) and ``seeded.tsv`` (read name, target: the
sample's pairs on them); and ``query.json``.
"""

import argparse
import json
import math
import mmap
import os
import sys
import time
from pathlib import Path
from typing import Final

import build
import numpy as np
import polars as pl
import query
import recruit
import tier2
import verify
from numpy.lib.format import open_memmap
from numpy.typing import NDArray

from kmer_functional_profiler import _core

U64 = np.uint64
PART_BITS: Final = 6  # 64 key partitions by code range: the merge sorts one at a time
PER_BUCKET: Final = 32  # postings per bucket aimed at (release: capped at 2^27 buckets)
MAX_BUCKET_BITS: Final = 27
MIX: Final = ((31, 0x94D049BB133111EB), (27, 0xBF58476D1CE4E5B9), (30, 1))


def _unshift(z: NDArray[np.uint64], s: int) -> NDArray[np.uint64]:
    """``x ^ (x >> s)``'s inverse."""
    x = z.copy()
    for _ in range(64 // s + 1):
        x = z ^ (x >> U64(s))
    return x


def codes(h: NDArray[np.uint64], k: int) -> NDArray[np.uint64]:
    """Each protein k-mer hash's base-20 rank: ``kfp::hash_kmer`` (splitmix64's finalizer, a
    bijection) undone to the packed k-mer (5 bits a residue), then read in base 20. Exact and
    dense: 39 bits at k 9."""
    z = np.asarray(h, dtype=U64)
    with np.errstate(over="ignore"):
        for s, mul in MIX:
            z = _unshift(z, s) * U64(pow(mul, -1, 1 << 64))
        packed = z - U64(0x9E3779B97F4A7C15)
        out = np.zeros(len(z), dtype=U64)
        for i in range(k):  # the first residue is the most significant
            out = out * U64(20) + ((packed >> U64(5 * (k - 1 - i))) & U64(31))
    return out


# a unit's representative, its first record: MGnify's cluster_rep, else (a slice without it) the
# longest member, as MMseqs picks centres (ties: the least id)
REP_ORDER: Final = (
    pl.col("cluster_rep"), pl.col("protein_id") != pl.col("cluster_rep"),
    pl.col("sequence").str.len_bytes(), pl.col("protein_id"),
)  # fmt: skip
REP_DESCENDING: Final = [False, False, True, False]


def unit_reps(members: pl.DataFrame) -> pl.DataFrame:
    """(protein_id, rep): each protein's unit representative, as ``shard`` orders them."""
    return (
        members.sort(*REP_ORDER, descending=REP_DESCENDING)
        .with_columns(rep=pl.col("protein_id").first().over("cluster_rep"))
        .select("protein_id", "rep")
    )


def code_bits(k: int) -> int:
    return math.ceil(k * math.log2(20))


# ---------------------------------------------------------------- build


def seqs(a: argparse.Namespace) -> None:
    """One shard of the sequences, built once for every index size: its units' proteins as
    FASTA, representative first, each unit's byte offset and its members."""
    out, s = Path(a.seqs), a.shard
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    m = (
        pl.scan_parquet(a.members)
        .select("protein_id", "cluster_rep", "sequence")
        .filter(pl.col("cluster_rep").hash(seed=1) % a.shards == s)
        .sort(*REP_ORDER, descending=REP_DESCENDING)
        .collect()
    )
    rep = m["cluster_rep"].to_numpy()
    first = np.flatnonzero(np.r_[True, rep[1:] != rep[:-1]]) if m.height else np.zeros(0, int)
    rec = (">" + m["protein_id"].cast(pl.String) + "\n" + m["sequence"] + "\n").to_list()
    start = np.r_[0, np.cumsum([len(r) for r in rec], dtype=np.int64)]
    with open(out / f"{s:04d}.fa", "w") as f:
        f.writelines(rec)
    np.save(out / f"{s:04d}.units.npy", start[np.r_[first, len(rec)]].astype(U64))
    np.save(out / f"{s:04d}.members.npy", np.diff(np.r_[first, m.height]).astype(np.uint32))
    got = {"shard": s, "shards": a.shards, "proteins": m.height, "units": len(first),
           "seqs_s": time.perf_counter() - t0}  # fmt: skip
    (out / f"{s:04d}.json").write_text(json.dumps(got) + "\n")
    print(json.dumps(got))


def shard_reps(d: Path, s: int) -> tuple[list[bytes], NDArray[np.int64]]:
    """A sequence shard's representatives (each unit's first record) and members per unit."""
    off = np.load(d / f"{s:04d}.units.npy").astype(np.int64)
    members = np.load(d / f"{s:04d}.members.npy").astype(np.int64)
    fa = (d / f"{s:04d}.fa").read_bytes()
    reps = []
    for lo in off[:-1].tolist():
        a = fa.index(b"\n", lo)
        reps.append(fa[a + 1 : fa.index(b"\n", a + 1)])
    return reps, members


def shard(a: argparse.Namespace) -> None:
    """One shard's keys at this index's size, by code partition, with shard-local units: the
    representatives read from the shared sequences (``seqs``)."""
    out, s = Path(a.index), a.shard
    t0 = time.perf_counter()
    reps, members = shard_reps(Path(a.seqs), s)
    got = {"shard": s, "units": len(reps), "keys": 0}
    (out / "shards").mkdir(parents=True, exist_ok=True)
    if not reps:  # a shard with no units (small references)
        (out / f"shards/{s:04d}.json").write_text(json.dumps(got) + "\n")
        return
    first = np.cumsum(members) - members  # each unit's representative: its first protein
    # keys: each representative's window minimizers, the rarest among the shard's per window
    info = build.kmer_info(reps, a.k, build.residue_log_freq(reps))
    vc = info.group_by("hash").len().sort("hash")
    rank = (vc["hash"].to_numpy(), vc["len"].to_numpy().astype(U64))
    unit, h = build.window_kmers(reps, a.k, a.window, 0, "protein", rank)
    mean = float(np.mean([max(len(r) - a.k + 1, 1) for r in reps])) / a.rate
    group = np.repeat(np.arange(len(reps)), members)  # each protein's unit
    pos, h = build.allocate(group, first[unit], h, f"sqrt:{mean}")
    keys = pl.DataFrame({"seq": group[pos], "hash": h}).join(
        info.select("seq", "hash", "bits0"), on=["seq", "hash"], how="left"
    )
    c = codes(keys["hash"].to_numpy(), a.k)
    keys = pl.DataFrame({
        "code": c, "unit": keys["seq"].cast(pl.UInt32),
        "bits0": keys["bits0"].fill_null(0.0).cast(pl.Float32),
        "part": (c >> U64(code_bits(a.k) - PART_BITS)).astype(np.int64),
    })  # fmt: skip
    for (p,), g in keys.partition_by("part", as_dict=True).items():
        d = out / f"parts/{p:02d}"
        d.mkdir(parents=True, exist_ok=True)
        g.drop("part").write_parquet(d / f"{s:04d}.parquet")
    got |= {"keys": keys.height, "mean": mean, "shard_s": time.perf_counter() - t0}
    (out / f"shards/{s:04d}.json").write_text(json.dumps(got) + "\n")
    print(json.dumps(got))


def merge(a: argparse.Namespace) -> None:
    """The shards' keys into one sorted posting table; each key's bits from its units."""
    out, sq = Path(a.index), Path(a.seqs)
    t0 = time.perf_counter()
    infos = sorted((json.loads(p.read_text()) for p in sq.glob("*.json")), key=lambda x: x["shard"])
    built = {json.loads(p.read_text())["shard"] for p in (out / "shards").glob("*.json")}
    if built != {x["shard"] for x in infos}:
        sys.exit(f"keys built for shards {sorted(built)}, sequences for {len(infos)}")
    base = np.r_[0, np.cumsum([x["units"] for x in infos])].astype(np.int64)
    np.save(out / "shards.npy", base)
    files = sorted((out / "parts").glob("*/*.parquet"))
    total = sum(pl.scan_parquet(f).select(pl.len()).collect().item() for f in files)
    b_bits = int(np.clip(round(math.log2(max(total / PER_BUCKET, 2))), 1, MAX_BUCKET_BITS))
    s_bits = code_bits(a.k) - b_bits
    keys = out / "keys"
    keys.mkdir(exist_ok=True)
    suffix = open_memmap(keys / "suffix.npy", "w+", np.uint16 if s_bits <= 16 else np.uint32,
                         (total,))  # fmt: skip
    unit = open_memmap(keys / "unit.npy", "w+", np.uint32, (total,))
    bits = open_memmap(keys / "bits.npy", "w+", np.uint8, (total,))
    count = np.zeros(2**b_bits, dtype=np.int64)
    at = 0
    for d in sorted((out / "parts").glob("*")):  # code ranges in order: the table stays sorted
        g = pl.concat([
            pl.read_parquet(f).with_columns(unit=pl.col("unit").cast(pl.Int64) + base[int(f.stem)])
            for f in sorted(d.glob("*.parquet"))
        ]).sort("code", "unit")  # fmt: skip
        g = g.with_columns(bits=pl.col("bits0") - pl.len().over("code").cast(pl.Float64).log(2))
        c, n = g["code"].to_numpy(), g.height
        suffix[at : at + n] = c & U64(2**s_bits - 1)
        unit[at : at + n] = g["unit"].to_numpy()
        bits[at : at + n] = query.encode_bits(g["bits"].to_numpy())
        count += np.bincount((c >> U64(s_bits)).astype(np.int64), minlength=len(count))
        at += n
    np.save(keys / "bucket.npy", np.r_[0, np.cumsum(count)].astype(U64))
    for x in (suffix, unit, bits):
        x.flush()
    meta = {"k": a.k, "window": a.window, "rate": a.rate, "bucket_bits": b_bits,
            "suffix_bits": s_bits, "postings": total, "units": int(base[-1]),
            "proteins": sum(x["proteins"] for x in infos), "shards": len(infos),
            "seqs": os.path.relpath(sq.resolve(), out.resolve()),
            "bytes": {f.name: f.stat().st_size for f in keys.glob("*.npy")},
            "merge_s": time.perf_counter() - t0}  # fmt: skip
    (out / "meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(json.dumps(meta))


# ---------------------------------------------------------------- query


class Compact:
    """The index, memory-mapped."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.meta = json.loads((self.path / "meta.json").read_text())
        self.k, self.s_bits = int(self.meta["k"]), int(self.meta["suffix_bits"])
        load = lambda n: np.load(self.path / f"keys/{n}.npy", mmap_mode="r")  # noqa: E731
        self.bucket, self.suffix, self.unit, self.bits = map(
            load, ("bucket", "suffix", "unit", "bits")
        )
        self.base = np.load(self.path / "shards.npy")
        self.seqs = self.path / self.meta["seqs"]  # the shared sequences (``seqs``)
        self._off: dict[int, NDArray[np.uint64]] = {}
        self._fa: dict[int, mmap.mmap] = {}

    def _search(self, s: NDArray[np.int64], lo: NDArray[np.int64], hi: NDArray[np.int64],
                right: bool) -> NDArray[np.int64]:  # fmt: skip
        """Per query, the first row of ``lo..hi`` whose suffix is >= s (> s if ``right``)."""
        lo, hi = lo.copy(), hi.copy()
        last = max(len(self.suffix) - 1, 0)
        while True:
            act = np.flatnonzero(lo < hi)
            if not len(act):
                return lo
            mid = (lo[act] + hi[act]) // 2
            v = np.asarray(self.suffix[np.minimum(mid, last)]).astype(np.int64)
            go = v <= s[act] if right else v < s[act]
            lo[act[go]] = mid[go] + 1
            hi[act[~go]] = mid[~go]

    def lookup(self, c: NDArray[np.uint64]) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
        """(query index, posting row) for each posting of each code."""
        b = (c >> U64(self.s_bits)).astype(np.int64)
        s = (c & U64(2**self.s_bits - 1)).astype(np.int64)
        lo, hi = self.bucket[b].astype(np.int64), self.bucket[b + 1].astype(np.int64)
        first = self._search(s, lo, hi, False)
        n = self._search(s, first, hi, True) - first
        q = np.repeat(np.arange(len(c)), n)
        return q, np.repeat(first - np.cumsum(n) + n, n) + np.arange(n.sum())

    def reps(self, units: NDArray[np.int64]) -> dict[int, tuple[int, bytes]]:
        """Each unit's representative: (protein id, sequence), its FASTA's first record."""
        out = {}
        sh = np.searchsorted(self.base, units, "right") - 1
        for s in np.unique(sh).tolist():
            if s not in self._off:
                self._off[s] = np.load(self.seqs / f"{s:04d}.units.npy")
                with open(self.seqs / f"{s:04d}.fa", "rb") as f:
                    self._fa[s] = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            off, fa = self._off[s], self._fa[s]
            for u in units[sh == s].tolist():
                lo, hi = int(off[u - self.base[s]]), int(off[u - self.base[s] + 1])
                rec = fa[lo : min(hi, lo + 65_536)]
                a = rec.find(b"\n")
                b = rec.find(b"\n", a + 1)
                if b < 0:  # a representative longer than the window
                    rec = fa[lo:hi]
                    b = rec.find(b"\n", a + 1)
                out[u] = (int(rec[1:a]), rec[a + 1 : b])
        return out


def tier1(ix: Compact, r1: str, r2: str) -> tuple[pl.DataFrame, float]:
    """(read, unit, score): each read pair's evidence on each unit its k-mers hit (bits summed
    less log2 windows per read pair for each hit beyond the first); and that charge."""
    parts, windows = [], 0
    it = _core.FastxHits(r1, r2, k=ix.k, frames="all")
    for b in it:
        c = codes(b["hash"], ix.k)
        windows += len(c)
        uc, inv = np.unique(c, return_inverse=True)
        q, _ = ix.lookup(uc)
        keyed = np.zeros(len(uc), dtype=bool)
        keyed[q] = True
        sel = keyed[inv]
        parts.append(pl.DataFrame({"read": b["read"][sel].astype(np.int64), "code": c[sel]}))
    log2_w = float(np.log2(max(windows / max(it.n_reads, 1), 2)))
    pairs = pl.concat(parts).unique() if parts else pl.DataFrame(
        schema={"read": pl.Int64, "code": pl.UInt64})  # fmt: skip
    uc = np.unique(pairs["code"].to_numpy())
    q, rows = ix.lookup(uc)
    post = pl.DataFrame({
        "code": uc[q], "unit": np.asarray(ix.unit[rows]).astype(np.int64),
        "bits": np.asarray(ix.bits[rows]).astype(np.float64) * query.BITS_STEP,
    })  # fmt: skip
    ev = (
        pairs.join(post, on="code")
        .group_by("read", "unit")
        .agg(score=pl.col("bits").sum() - (pl.len() - 1) * log2_w)
    )
    return ev, log2_w


def confirm(ix: Compact, reads: list[str], seeded: pl.DataFrame,
            a: argparse.Namespace) -> pl.DataFrame:  # fmt: skip
    """(unit, score): each unit's best seeded read pair against its representative, a read
    batch at a time (only that batch's representatives in memory: at release most seeded
    units are chance hits, ~8 per read pair)."""
    out = []
    for lo, m1, m2 in recruit.batches(*reads, a.read_batch):
        mine = seeded.filter(pl.col("read").is_between(lo, lo + len(m1) - 1))
        r = np.unique(mine["read"].to_numpy()).tolist()
        fr = _core.pair_frames([m1[i - lo] for i in r], [m2[i - lo] for i in r])
        seqs = {
            u: verify.mask(x) for u, (_, x) in ix.reps(mine["unit"].unique().to_numpy()).items()
        }
        got = recruit.aligned(mine.rename({"unit": "target"}), dict(zip(r, fr, strict=True)),
                              seqs, a.batch_pairs).rename({"target": "unit"})  # fmt: skip
        out.append(got.group_by("unit").agg(pl.col("score").max()))
    return pl.concat(out).group_by("unit").agg(pl.col("score").max())


def run_query(a: argparse.Namespace) -> None:
    ix, out = Compact(a.index), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    info: dict[str, float | int] = {}
    reads = {"sample": a.reads, "decoy": query.reversed_reads(a.reads, out)}
    seeded, best = {}, {}
    t = time.perf_counter()
    for name, (r1, r2) in reads.items():
        ev, info[f"{name}_log2_windows"] = tier1(ix, r1, r2)
        seeded[name] = ev.filter(pl.col("score") >= a.floor)
        info[f"{name}_tier1_pairs"] = seeded[name].height
    info["tier1_s"] = time.perf_counter() - t
    # confirmation: each seeded pair against its unit's representative
    t = time.perf_counter()
    cand = pl.concat([x.select("unit") for x in seeded.values()])["unit"].n_unique()
    for name, r in reads.items():
        best[name] = confirm(ix, r, seeded[name], a)
    bar, _ = query.decoy_bar(best["sample"], best["decoy"], a.share, scale=a.chance_scale)
    confirmed = np.sort(best["sample"].filter(pl.col("score") >= bar)["unit"].to_numpy())
    info |= {"candidate_units": cand, "confirm_bar": bar, "confirmed_units": len(confirmed),
             "confirm_s": time.perf_counter() - t}  # fmt: skip
    # the confirmed units for recruitment (recruit.py, or any other recruiting tool): their
    # representatives, and the sample's pairs on them
    rep = ix.reps(confirmed)
    pid = {u: rep[u][0] for u in confirmed.tolist()}
    with open(out / "confirmed.fa", "wb") as f:
        f.writelines(b">%d\n%s\n" % (pid[u], rep[u][1]) for u in confirmed.tolist())
    names = recruit.read_names(a.reads[0])
    on = seeded["sample"].filter(pl.col("unit").is_in(pl.Series(confirmed).implode()))
    pl.DataFrame({
        "read": pl.Series(names).gather(on["read"]), "target": on["unit"].replace_strict(pid)
    }).write_csv(out / "seeded.tsv", separator="\t")  # fmt: skip
    best["sample"].filter(pl.col("score") >= bar).with_columns(
        protein_id=pl.col("unit").replace_strict(pid)
    ).sort("unit").write_csv(out / "units.tsv", separator="\t", float_precision=2)
    for f in out.glob("decoy_*.fq.gz"):
        f.unlink()
    info["peak_rss_gb"] = tier2.peak_rss_gb()
    (out / "query.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info))


def main() -> None:
    p = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    sq = sub.add_parser("seqs", help="one shard of the sequences, shared by every index size")
    sq.add_argument("--members", required=True, help="parquet: protein_id, cluster_rep, sequence")
    sq.add_argument("--shards", type=int, default=1)
    s = sub.add_parser("shard", help="one shard's keys at this index's size")
    m = sub.add_parser("merge", help="the shards' keys into the posting table")
    for x in (sq, s):
        x.add_argument("--shard", type=int, default=0)
    for x in (sq, s, m):
        x.add_argument("--seqs", required=True, help="the shared sequences' directory")
    for x in (s, m):
        x.add_argument("--index", required=True)
        x.add_argument("--k", type=int, default=9)
        x.add_argument("--window", type=int, default=4, help="window minimizers' window (k-mers)")
        x.add_argument("--rate", type=float, default=32, help="one key per this many rep residues")
    q = sub.add_parser("query", help="tier 1 and confirmation")
    q.add_argument("--index", required=True)
    q.add_argument("--reads", nargs=2, required=True)
    q.add_argument("--out", required=True)
    q.add_argument("--floor", type=float, default=30.0, help="tier-1 bits a read pair needs")
    q.add_argument("--share", type=float, default=0.01, help="decoy share at the bars")
    q.add_argument("--chance-scale", type=float, default=1.0, help="query.decoy_bar's scale")
    q.add_argument("--read-batch", type=int, default=200_000)
    q.add_argument("--batch-pairs", type=int, default=2_000_000)
    a = p.parse_args()
    {"seqs": seqs, "shard": shard, "merge": merge, "query": run_query}[a.cmd](a)


if __name__ == "__main__":
    main()
