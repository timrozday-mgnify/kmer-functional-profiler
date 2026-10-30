"""Tree index, phase 1: build one tier (docs/tree-index-plan.md, Implementation plan, phase 1,
with the survey's floor design and the proposed tiers).

The topology is the survey's (``survey.Tree``: the exact storage tree on per-sequence floored
sets at the topology rates). Leaves are numbered by their position in its leaf order, so every
node is a range of positions, ``lo..hi``. Postings are range postings, ``(hash, lo, hi, sign)``:
a leaf holds a k-mer iff the signs of the postings covering it sum to 1 (they sum to 0 or 1).

- **base** (kind 0): each leaf's representative's k-mers with hash <= t\\_base (and with
  ``--window`` W, its least hash in every window of W k-mers, any hash), diffed by Fitch over
  the whole tree: a gain is +1 over the node's range, a loss -1. With ``--ladder``
  k:max\\_len:w,..., leaves shorter than max\\_len also hold window minimizers of k-mers of that
  shorter k (hashes salted per k), for reads that overlap a short gene by too little to hold a
  whole k-mer of the index's k;
- **floor** (kind 1): per piece of two or more leaves (``survey.pieces`` at ``tau``), its
  n\\_min floor k-mers, +1 over the piece. Single-leaf pieces get no floor (*n\\_single* 0);
- **top-up** (kind 2): each leaf topped up to ``n_1`` k-mers of its own (base, floor k-mers it
  holds, then its lowest other hashes, beyond t\\_cap if need be), +1 over the leaf. Proteins
  with no k-mer up to t\\_cap get a position of their own after the tree's, holding their n\\_1
  lowest hashes.

Writes to ``--out``: ``hash.npy`` (u64), ``lo.npy``, ``hi.npy`` (u32), ``sign.npy`` (i8) and
``kind.npy`` (u8), sorted by (hash, lo); ``node_lo.npy``, ``node_hi.npy``, ``node_level.npy``
(internal node i of the merge tree, by positions); ``pos_size.npy`` and ``node_size.npy``
(|S\\_v|: k-mers held at each position and each internal node); ``pos_len.npy`` (each
position's representative's length); ``edge_a.npy``, ``edge_b.npy``, ``edge_w.npy`` (the spanning
tree's edges between positions, weighed in topology units, for the neighbour rule);
``pieces.parquet`` (lo, hi, floor,
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
import query
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
    extra_pos, extra_h = np.zeros(0, np.int64), np.zeros(0, U64)
    if a.window:
        reps = [tree.seqs[i] for i in tree.rep_by_pos]
        wp, wh = window_kmers(reps, a.k, a.window, a.window_min)
        rows = pl.DataFrame({"pos": tree.pos, "h": tree.h}).with_row_index()
        win = rows.join(pl.DataFrame({"pos": wp, "h": wh}), on=["pos", "h"], how="semi")
        keep[win["index"].to_numpy()] = True
        beyond = wh > U64(_core.max_hash(a.t_cap))  # not among the rows, which stop at t_cap
        extra_pos, extra_h = wp[beyond], wh[beyond]
    lengths = np.array([len(tree.seqs[i]) for i in tree.rep_by_pos])
    for k, max_len, w in query.ladder(a.ladder):  # short leaves: shorter k-mers too
        short = np.flatnonzero(lengths < max_len)
        sp, sh = window_kmers([tree.seqs[i] for i in tree.rep_by_pos[short]], k, w, 0)
        extra_pos = np.concatenate([extra_pos, short[sp]])
        extra_h = np.concatenate([extra_h, query.salted(sh, k)])
    pos = np.concatenate([tree.pos[keep], extra_pos])
    h = np.concatenate([tree.h[keep], extra_h])
    leaf = lay["order"][pos]
    o = np.lexsort((h, leaf))
    # every hash a leaf does not hold is known absent (no missing data above t_base)
    top = U64(ALL) if a.window or a.ladder else lo_base
    sets = p0.Sets(p0.csr(leaf[o], h[o], n), h[o], np.zeros(n, np.uint8), np.array([top]))
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
    spare = ~keep
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


def window_kmers(
    seqs: list[bytes], k: int, w: int, min_count: int = 0
) -> tuple[NDArray[np.int64], NDArray[np.uint64]]:
    """Distinct (sequence, hash): each sequence's least hash in every window of w_s consecutive
    k-mers (robust winnowing; the whole sequence when shorter), so every stretch of w_s + k - 1
    residues holds an indexed k-mer. w_s is ``w``, or with ``min_count`` M, min(w, n_k // M)
    (at least 1), so a sequence of n_k k-mers gets about M or more."""
    out_seq, out_h = [], []
    none = np.iinfo(np.uint64).max
    for start in range(0, len(seqs), 20_000):
        batch = seqs[start : start + 20_000]
        n_k = np.array([max(len(s) - k + 1, 0) for s in batch])
        # ponytail: each k-mer hashed as its own sequence; a positional kernel at scale
        o = _core.hash_proteins([s[i : i + k] for s in batch for i in range(len(s) - k + 1)], k)
        h = np.full(n_k.sum(), none, dtype=np.uint64)
        h[o["seq"]] = o["hash"]
        first = np.cumsum(n_k) - n_k
        w_s = np.clip(n_k // min_count, 1, w) if min_count else np.full(len(batch), w)
        for ww in np.unique(w_s[n_k > 0]).tolist():
            sel = np.flatnonzero((w_s == ww) & (n_k > 0))
            idx = np.concatenate([np.arange(first[i], first[i] + n_k[i]) for i in sel])
            owner, m = winnow(h[idx], n_k[sel], ww)
            out_seq.append(sel[owner] + start)
            out_h.append(m)
    seq = np.concatenate(out_seq) if out_seq else np.zeros(0, np.int64)
    hh = np.concatenate(out_h) if out_h else np.zeros(0, np.uint64)
    pairs = pl.DataFrame({"seq": seq, "h": hh}).filter(pl.col("h") != none).unique()
    return pairs["seq"].to_numpy().astype(np.int64), pairs["h"].to_numpy()


def winnow(
    h: NDArray[np.uint64], n_k: NDArray[np.int64], w: int
) -> tuple[NDArray[np.int64], NDArray[np.uint64]]:
    """(sequence, least hash) of every window of ``w`` within each sequence of ``h``
    (concatenated, ``n_k`` each); a sequence shorter than ``w`` is one window."""
    owner = np.repeat(np.arange(len(n_k)), n_k)
    out_o, out_m = [], []
    if len(h) >= w:
        m = np.lib.stride_tricks.sliding_window_view(h, w).min(axis=1)
        i = np.arange(len(m))
        ok = owner[i] == owner[i + w - 1]
        out_o.append(owner[i[ok]])
        out_m.append(m[ok])
    short = np.flatnonzero(n_k < w)
    if len(short):
        whole = np.minimum.reduceat(h, np.cumsum(n_k) - n_k)
        out_o.append(short)
        out_m.append(whole[short])
    return np.concatenate(out_o), np.concatenate(out_m)


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
        "--window", type=int, default=0,
        help="also base: each leaf's least hash in every window of this many k-mers (0: off)",
    )  # fmt: skip
    parser.add_argument(
        "--ladder", default="",
        help="k:max_len:window,...: leaves shorter than max_len also hold window minimizers of k",
    )  # fmt: skip
    parser.add_argument(
        "--window-min", type=int, default=0,
        help="windows shrink to n // this for a leaf of n k-mers, so short ones get more (0: off)",
    )  # fmt: skip
    parser.add_argument(
        "--leaves", choices=["proteins", "collapsed"], default="proteins",
        help="proteins: a leaf per distinct dense set; collapsed: phase 0's (survey.Tree)",
    )  # fmt: skip
    parser.add_argument("--topology-tree", choices=["spanning", "merge"], default="spanning")
    parser.add_argument("--topology-k", type=int, default=8)
    parser.add_argument("--topology-t", type=float, default=1.0)
    parser.add_argument("--topology-max-holders", type=int, default=0)
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
    pos_of_leaf = np.empty(n, dtype=np.int64)
    pos_of_leaf[lay["order"]] = np.arange(n)
    f = tree.forest
    np.save(out / "edge_a.npy", pos_of_leaf[f.mst_a].astype(np.uint32))
    np.save(out / "edge_b.npy", pos_of_leaf[f.mst_b].astype(np.uint32))
    np.save(out / "edge_w.npy", f.mst_w.astype(np.uint32))
    lengths = np.array([len(s) for s in tree.seqs], dtype=np.uint32)
    leafless = np.flatnonzero((tree.leaf_of < 0) & (tree.n_valid > 0)) if a.n_1 else []
    np.save(out / "pos_len.npy", np.concatenate([lengths[tree.rep_by_pos], lengths[leafless]]))
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
            "kruskal_s": tree.kruskal_s,
            "spanning_tree_s": tree.spanning_tree_s,
            "build_s": time.perf_counter() - t0,
            "peak_rss_bytes": rss if sys.platform == "darwin" else rss * 1024,
        },
        "files": files,
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(json.dumps(meta["stats"]))


if __name__ == "__main__":
    main()
