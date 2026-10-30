"""Tree index, phase 1: build one tier (docs/tree-index-plan.md, Implementation plan, phase 1,
with the survey's floor design and the proposed tiers).

The topology is the survey's (``survey.Tree``: the exact storage tree on per-sequence floored
sets at the topology rates). Leaves are numbered by their position in its leaf order, so every
node is a range of positions, ``lo..hi``. Postings are range postings, ``(hash, lo, hi, sign)``:
a leaf holds a k-mer iff the signs of the postings covering it sum to 1 (they sum to 0 or 1).

- **base** (kind 0): each leaf's representative's k-mers with hash <= t\\_base, diffed by Fitch
  over the whole tree: a gain is +1 over the node's range, a loss -1;
- **floor** (kind 1): per piece of two or more leaves (``survey.pieces`` at ``tau``), its
  n\\_min floor k-mers, +1 over the piece. Single-leaf pieces get no floor (*n\\_single* 0);
- **top-up** (kind 2): each leaf topped up to ``n_1`` k-mers of its own (base, floor k-mers it
  holds, then its lowest other hashes, beyond t\\_cap if need be), +1 over the leaf. Proteins
  with no k-mer up to t\\_cap get a position of their own after the tree's, holding their n\\_1
  lowest hashes.

Writes to ``--out``: ``hash.npy`` (u64), ``lo.npy``, ``hi.npy`` (u32), ``sign.npy`` (i8) and
``kind.npy`` (u8), sorted by (hash, lo); ``node_lo.npy``, ``node_hi.npy``, ``node_level.npy``
(internal node i of the merge tree, by positions); ``pos_size.npy`` and ``node_size.npy``
(|S\\_v|: k-mers held at each position and each internal node); ``pieces.parquet`` (lo, hi, floor,
coherence); ``proteins.parquet`` (protein\\_id, cluster\\_rep, pos; -1 for no valid k-mer); and
``meta.json``.
"""

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import phase0 as p0
import polars as pl
import survey
from numpy.typing import NDArray

from kmer_functional_profiler import _core

U64 = np.uint64
ALL = 2**64 - 1
BASE, FLOOR, TOPUP = 0, 1, 2


def postings(tree: survey.Tree, a: argparse.Namespace) -> tuple[pl.DataFrame, pl.DataFrame, int]:
    """Range postings of every kind, the pieces, and positions used (tree's plus leafless)."""
    f, lay, n = tree.forest, tree.lay, tree.forest.n
    lo_base = U64(_core.max_hash(a.t_base))
    # base band: Fitch's changes over the whole tree
    keep = tree.h <= lo_base
    pos, h = tree.pos[keep], tree.h[keep]
    leaf = lay["order"][pos]
    o = np.lexsort((h, leaf))
    sets = p0.Sets(p0.csr(leaf[o], h[o], n), h[o], np.zeros(n, np.uint8), np.array([lo_base]))
    ch = _core.tree_schemes(
        n, f.left, f.right, f.tops, sets.offsets.astype(U64), sets.hashes, sets.level,
        sets.level_hash, emit=True,
    )  # fmt: skip
    node = ch["change_node"]
    parts = [
        pl.DataFrame(
            {
                "hash": ch["change_hash"],
                "lo": lay["lo"][node],
                "hi": lay["hi"][node],
                "sign": np.where(ch["change_gain"], 1, -1).astype(np.int8),
                "kind": np.uint8(BASE),
            }
        )
    ]
    own = np.bincount(pos, minlength=n)  # k-mers each position holds of its own
    # floors per piece; rows above t_base that are not floor k-mers can top a leaf up
    got, _ = survey.pieces(tree, a.t_base, a.n_min, a.tau)
    spare = tree.h > lo_base
    rows = []
    for lo, hi, floor, coh in got:
        rows.append((lo, hi, len(floor) if hi - lo > 1 else 0, coh))
        if hi - lo == 1:
            continue
        parts.append(
            pl.DataFrame(
                {
                    "hash": floor,
                    "lo": np.uint32(lo),
                    "hi": np.uint32(hi),
                    "sign": np.int8(1),
                    "kind": np.uint8(FLOOR),
                }
            )  # fmt: skip
        )
        r0, r1 = tree.row[lo], tree.row[hi]
        hit = spare[r0:r1] & np.isin(tree.h[r0:r1], floor)
        spare[r0:r1] &= ~hit
        own[lo:hi] += np.bincount(tree.pos[r0:r1][hit] - lo, minlength=hi - lo)
    pieces = pl.DataFrame(rows, schema=["lo", "hi", "floor", "coherence"], orient="row")
    # top-ups: each position's lowest spare hashes (rows are sorted by hash within a position)
    need = np.maximum(a.n_1 - own, 0)
    sp, sh = tree.pos[spare], tree.h[spare]
    rank = np.arange(len(sp)) - np.searchsorted(sp, sp)
    take = rank < need[sp]
    parts.append(topups(sh[take], sp[take]))
    short = need - np.bincount(sp[take], minlength=n)
    # beyond t_cap where up to t_cap is not enough
    lack = np.flatnonzero(short > 0)
    beyond = _core.max_hash(a.t_cap) + 1
    if len(lack):
        reps = [tree.seqs[i] for i in tree.rep_by_pos[lack]]
        seq, hh = p0.hashed(reps, a.k, ALL, min_hash=beyond)
        rank = np.arange(len(seq)) - np.searchsorted(seq, seq)
        take = rank < short[lack][seq]
        parts.append(topups(hh[take], lack[seq[take]]))
    # proteins without a leaf: positions of their own after the tree's
    leafless = np.flatnonzero((tree.leaf_of < 0) & (tree.n_valid > 0)) if a.n_1 else []
    if len(leafless):
        seq, hh = p0.hashed([tree.seqs[i] for i in leafless], a.k, ALL)
        rank = np.arange(len(seq)) - np.searchsorted(seq, seq)
        take = rank < a.n_1
        parts.append(topups(hh[take], n + seq[take]))
    table = pl.concat(parts).sort("hash", "lo", "hi")
    return table, pieces, n + len(leafless)


def topups(h: NDArray[np.uint64], pos: NDArray[np.int64]) -> pl.DataFrame:
    return pl.DataFrame(
        {"hash": h, "lo": pos.astype(np.uint32), "hi": (pos + 1).astype(np.uint32),
         "sign": np.int8(1), "kind": np.uint8(TOPUP)}
    )  # fmt: skip


def protein_positions(tree: survey.Tree, n_1: int) -> NDArray[np.int64]:
    n = tree.forest.n
    pos_of_leaf = np.empty(n, dtype=np.int64)
    pos_of_leaf[tree.lay["order"]] = np.arange(n)
    out = np.full(tree.n_proteins, -1, dtype=np.int64)
    has = tree.leaf_of >= 0
    out[has] = pos_of_leaf[tree.leaf_of[has]]
    if n_1:
        leafless = np.flatnonzero(~has & (tree.n_valid > 0))
        out[leafless] = n + np.arange(len(leafless))
    return out


def holdings(index: Path) -> pl.DataFrame:
    """Every (hash, pos) the index holds, by expanding its ranges (brute force, for checks).
    Raises if a cover sums to anything but 0 or 1."""
    t = pl.DataFrame({x: np.load(index / f"{x}.npy") for x in ("hash", "lo", "hi", "sign")})
    cover = (
        t.with_columns(pos=pl.int_ranges("lo", "hi"))
        .explode("pos")
        .group_by("hash", "pos")
        .agg(pl.col("sign").cast(pl.Int32).sum())
    )
    bad = cover.filter(~pl.col("sign").is_in([0, 1]))
    if bad.height:
        raise ValueError(f"{bad.height} covers sum outside 0..1, e.g. {bad.row(0)}")
    return cover.filter(pl.col("sign") == 1).select("hash", "pos").sort("pos", "hash")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--members", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--k", type=int, default=11)
    parser.add_argument("--t-cap", type=float, default=0.2)
    parser.add_argument("--levels", type=int, default=8)
    parser.add_argument("--topology-t-base", type=float, default=0.02)
    parser.add_argument("--topology-n-min", type=int, default=16)
    parser.add_argument("--t-base", type=float, default=0.02)
    parser.add_argument("--n-min", type=int, default=16)
    parser.add_argument("--tau", type=float, default=0.3)
    parser.add_argument("--n-1", type=int, default=2)
    parser.add_argument(
        "--leaves", choices=["proteins", "collapsed"], default="proteins",
        help="proteins: a leaf per distinct protein set; collapsed: phase 0's (survey.Tree)",
    )  # fmt: skip
    a = parser.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    members = p0.load_members(a.members)
    tree = survey.Tree(members, a)
    topology_s = time.perf_counter() - t0
    table, pieces, n_pos = postings(tree, a)
    for col, dtype in (("hash", U64), ("lo", np.uint32), ("hi", np.uint32), ("sign", np.int8),
                       ("kind", np.uint8)):  # fmt: skip
        np.save(out / f"{col}.npy", table[col].to_numpy().astype(dtype))
    lay, n = tree.lay, tree.forest.n
    np.save(out / "node_lo.npy", lay["lo"][n:])
    np.save(out / "node_hi.npy", lay["hi"][n:])
    np.save(out / "node_level.npy", tree.forest.level)
    # |S_v|: k-mers held at each position and each internal node (Fitch's states)
    rng = (
        table["lo"].to_numpy(),
        table["hi"].to_numpy(),
        table["sign"].to_numpy().astype(np.int64),
    )
    pos = np.arange(n_pos, dtype=np.uint32)
    np.save(out / "pos_size.npy", _core.contained_sums(*rng, pos, pos + 1).astype(np.uint32))
    size = _core.contained_sums(*rng, lay["lo"][n:], lay["hi"][n:])
    np.save(out / "node_size.npy", size.astype(np.uint32))
    pieces.write_parquet(out / "pieces.parquet")
    members.select("protein_id", "cluster_rep").with_columns(
        pos=protein_positions(tree, a.n_1)
    ).write_parquet(out / "proteins.parquet")
    kinds = np.bincount(table["kind"].to_numpy(), minlength=3)
    files = {p.name: p.stat().st_size for p in sorted(out.iterdir()) if p.name != "meta.json"}
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    meta = {
        "format": "tree-v0",
        "params": {k: v for k, v in vars(a).items() if k not in ("members", "out")},
        "stats": {
            "proteins": tree.n_proteins,
            "units": int(members["cluster_rep"].n_unique()),
            "positions": n_pos,
            "tree_leaves": n,
            "internal_nodes": len(tree.forest.left),
            "components": len(tree.forest.tops),
            "pieces": pieces.height,
            "floored_pieces": int((pieces["floor"] > 0).sum()),
            "postings": table.height,
            "base_postings": int(kinds[BASE]),
            "floor_postings": int(kinds[FLOOR]),
            "topup_postings": int(kinds[TOPUP]),
            "keys": table["hash"].n_unique(),
            "bytes": sum(files.values()),
            "bytes_per_posting": sum(files.values()) / max(table.height, 1),
            "topology_s": topology_s,
            "build_s": time.perf_counter() - t0,
            "peak_rss_bytes": rss if sys.platform == "darwin" else rss * 1024,
        },
        "files": files,
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(json.dumps(meta["stats"]))


if __name__ == "__main__":
    main()
