"""Tree index: rough alignment of candidates, a trial (docs/tree-index-plan.md, Progress log:
gather, per-layer bars and layer representatives).

Each (read pair, leaf) with evidence (``query.evidence`` per read pair, at least ``--floor``
bits) is scored by the best ungapped BLOSUM62 diagonal of the read pair's translated frames
that hold a key (of both mates' six, joined by stops) against the leaf's representative protein (its
longest member, as ``build.py``'s), in bits (Karlin-Altschul, ungapped lambda and K). The
reads reversed are scored the same way, for target-decoy bars at each ``--shares``:

- ``verify_<share>``: leaves whose best read pair scores at least the bar;
- ``vgather_<share>``: ``query.gather`` over the read pairs scoring at least that bar: the
  fewest leaves that explain every verified read pair (fewer targets for tier 2).

Writes ``cand_verify_<share>.npy`` and ``cand_vgather_<share>.npy`` into ``--query`` (for
``recall.py score``), ``aligned.parquet`` (pos, read, score: the sample's pairs) and
``verify.json`` (pairs, time, pairs per second, bars).
"""

import argparse
import gzip
import json
import math
import time
from pathlib import Path

import numpy as np
import polars as pl
import query
from numpy.typing import NDArray

from kmer_functional_profiler import _core

LAMBDA, K = 0.3176, 0.134  # BLOSUM62 ungapped (BLAST)


def blosum62() -> NDArray[np.int8]:
    """BLOSUM62 by byte pair (``a << 8 | b``); any pair with a byte outside it, stops
    included, -20, so a segment does not run across frames."""
    from Bio.Align import substitution_matrices

    m = substitution_matrices.load("BLOSUM62")
    out = np.full(1 << 16, -20, dtype=np.int8)
    for x in m.alphabet:
        for y in m.alphabet:
            if x != "*" and y != "*":
                out[ord(x) << 8 | ord(y)] = int(m[x][y])
    return out


def read_frames(index: query.Index, r1: str, r2: str) -> list[bytes]:
    """Per read pair (in file order, as ``FastxHits`` numbers them), the translated frames
    (mate and frame) holding a key at any of the index's k, joined by stops."""
    keyed = []
    todo = [(index.k, "protein", 0)] + [(k, al, c) for k, _, al, c in query.steps(index)]
    for k, alpha, code in todo:
        for b in _core.FastxHits(r1, r2, k=k, frames="all", alphabet=alpha):
            h = query.salted(b["hash"], code) if code else b["hash"]
            hit = np.isin(h, index.keyed(h))
            keyed.append(
                b["read"][hit].astype(np.int64) * 12 + b["mate"][hit] * 6 + b["frame"][hit]
            )
    want = set(np.unique(np.concatenate(keyed)).tolist())
    out = []
    with gzip.open(r1, "rb") as f, gzip.open(r2, "rb") as g:
        for i, (a, b) in enumerate(zip(f, g, strict=True)):
            if i % 4 == 1:
                r = i // 4
                frames = _core.translate_frames(a.strip()) + _core.translate_frames(b.strip())
                out.append(b"*".join(x for j, x in enumerate(frames) if r * 12 + j in want))
    return out


def rep_seqs(index: Path, members: str) -> dict[int, bytes]:
    """Each position's representative: its longest protein (ties: the first)."""
    prot = pl.read_parquet(index / "proteins.parquet").filter(pl.col("pos") >= 0)
    seqs = pl.read_parquet(members, columns=["protein_id", "sequence"])
    got = (
        prot.join(seqs, on="protein_id")
        .with_columns(n=pl.col("sequence").str.len_bytes())
        .sort("pos", "n", descending=[False, True])
        .unique("pos", keep="first")
    )
    return dict(zip(got["pos"].to_list(), (s.encode() for s in got["sequence"]), strict=True))


def scored(
    ev: pl.DataFrame, frames: list[bytes], reps: dict[int, bytes], matrix: NDArray[np.int8]
) -> tuple[pl.DataFrame, float]:
    """(pos, read, score): each pair's alignment bits; and the seconds the alignment took."""
    reads, pos = ev["read"].to_list(), ev["pos"].to_list()
    t0 = time.perf_counter()
    s = _core.align_ungapped([frames[r] for r in reads], [reps[p] for p in pos], matrix)
    took = time.perf_counter() - t0
    bits = (LAMBDA * s.astype(np.float64) - math.log(K)) / math.log(2)
    return ev.select("pos", "read").with_columns(score=bits), took


def main() -> None:
    p = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    p.add_argument("--index", required=True)
    p.add_argument("--members", required=True, help="the index's members (sequences)")
    p.add_argument("--reads", nargs=2, required=True)
    p.add_argument("--query", required=True, help="query.py --out: candidates are added here")
    p.add_argument("--floor", type=float, default=8.0, help="read-pair evidence to align")
    p.add_argument("--shares", default="0.01,0.05,0.1")
    a = p.parse_args()
    index, qdir = query.Index(a.index), Path(a.query)
    matrix, reps = blosum62(), rep_seqs(Path(a.index), a.members)
    info: dict[str, float | int] = {}
    got = {}
    for name, reads in (("sample", a.reads), ("decoy", query.reversed_reads(a.reads, qdir))):
        t0 = time.perf_counter()
        pairs, windows, n_reads = query.read_pairs(index, *reads, "all")
        log2_w = float(np.log2(max(windows.get(index.k, 1) / max(n_reads, 1), 2)))
        ev = query.evidence_chunked(index, pairs, log2_w, True, per_read=True, floor=a.floor)
        ev = ev.filter(pl.col("pos").is_in(list(reps)))  # leafless positions: none here
        frames = read_frames(index, *reads)
        got[name], took = scored(ev, frames, reps, matrix)
        cells = sum(len(frames[r]) * len(reps[q]) for r, q in ev.select("read", "pos").iter_rows())
        info |= {f"{name}_pairs": ev.height, f"{name}_align_s": took,
                 f"{name}_pairs_per_s": ev.height / max(took, 1e-9),
                 f"{name}_cells_per_s": cells / max(took, 1e-9),
                 f"{name}_total_s": time.perf_counter() - t0}  # fmt: skip
    got["sample"].write_parquet(qdir / "aligned.parquet")
    best = {n: g.group_by("pos").agg(pl.col("score").max()) for n, g in got.items()}
    for s in map(float, a.shares.split(",")):
        bar, _ = query.decoy_bar(best["sample"], best["decoy"], s)
        np.save(
            qdir / f"cand_verify_{s}.npy",
            np.sort(best["sample"].filter(pl.col("score") >= bar)["pos"].to_numpy()),
        )
        # gather over the pairs at or above the bar: the fewest leaves explaining them
        g = query.gather(got["sample"].filter(pl.col("score") >= bar), 0.0)
        np.save(qdir / f"cand_vgather_{s}.npy", np.sort(g["pos"].to_numpy()))
        info |= {f"verify_{s}_bar": bar}
    for f in qdir.glob("decoy_*.fq.gz"):
        f.unlink()
    (qdir / "verify.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info))


if __name__ == "__main__":
    main()
