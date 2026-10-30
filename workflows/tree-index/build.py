"""Tree index, phase 1: build one tier (docs/tree-index-plan.md, Implementation plan, phase 1,
with the survey's floor design and the proposed tiers).

The topology is the survey's (``survey.Tree``: the exact storage tree on per-sequence floored
sets at the topology rates). Leaves are numbered by their position in its leaf order, so every
node is a range of positions, ``lo..hi``. Postings are range postings, ``(hash, lo, hi, sign)``:
a leaf holds a k-mer iff the signs of the postings covering it sum to 1 (they sum to 0 or 1).

- **base** (kind 0): each leaf's representative's k-mers with hash <= t\\_base (and with
  ``--window`` W, its least hash in every window of W k-mers, any hash), diffed by Fitch over
  the whole tree: a gain is +1 over the node's range, a loss -1. With ``--ladder``
  k:w,..., a leaf with fewer than ``--min-kmers`` informative k-mers (a priori rarer than one
  chance holder among the leaves: -log2 q0 >= log2 N) at the k before also holds window
  minimizers of that shorter k (hashes salted per k);
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
``keys.npy``, ``key_bits.npy`` and ``key_k.npy`` (each key's evidence per hit in bits,
-log2 q0 - log2 g, and its k); ``pieces.parquet`` (lo, hi, floor,
coherence); ``proteins.parquet`` (protein\\_id, cluster\\_rep, pos; -1 for no valid k-mer); and
``meta.json``.
"""

import argparse
import json
import resource
import sys
import time
from collections.abc import Callable
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


def postings(
    tree: survey.Tree, a: argparse.Namespace
) -> tuple[pl.DataFrame, pl.DataFrame, int, pl.DataFrame]:
    """Range postings of every kind, the pieces, positions used (tree's plus leafless), and
    each key's evidence (:func:`key_bits`)."""
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
    # bits of composition per k-mer (-log2 q0), and the ladder: a leaf takes the next k while it
    # has fewer than min_kmers informative k-mers (a priori rarer than one chance holder among
    # the leaves) at the k before
    n_ladder = 0  # ladder rows: the last of extra_*
    log_f = residue_log_freq(tree.seqs)
    reps = [tree.seqs[i] for i in tree.rep_by_pos]
    info = [kmer_info(tree.seqs, a.k, log_f).with_columns(k=pl.lit(a.k, pl.UInt8))]
    rare = np.log2(max(n, 2))
    need = np.flatnonzero(informative(kmer_info(reps, a.k, log_f), n, rare) < a.min_kmers)
    for k, w in query.ladder(a.ladder):
        if not len(need):
            break
        sub = [reps[i] for i in need]
        sp, sh = window_kmers(sub, k, w, a.window_min)
        extra_pos = np.concatenate([extra_pos, need[sp]])
        extra_h = np.concatenate([extra_h, query.salted(sh, k)])
        n_ladder += len(sp)
        got = kmer_info(sub, k, log_f)
        info.append(got.with_columns(hash=pl.Series(query.salted(got["hash"].to_numpy(), k)),
                                     k=pl.lit(k, pl.UInt8)))  # fmt: skip
        need = need[informative(got, len(need), rare) < a.min_kmers]
    pos = np.concatenate([tree.pos[keep], extra_pos])
    h = np.concatenate([tree.h[keep], extra_h])
    if a.detect_aai or a.unit_kmers:  # units: members share k-mers, which replace theirs
        if a.unit_kmers:  # dense intersections, n spread shared k-mers per unit
            dense = kmer_info(reps, a.k, log_f).sort("seq", "hash")
            d_row = np.searchsorted(dense["seq"].to_numpy(), np.arange(n + 1))
            d_h = dense["hash"].to_numpy()
            held = lambda v: d_h[d_row[v] : d_row[v + 1]]  # noqa: E731
            if a.unit_tau < 1:
                unit, shared = soft_units(tree, a.unit_kmers, a.unit_tau, held)
            else:
                unit, shared = detect_units(tree, a.unit_kmers, held)
            chosen = {u: spread_shared(reps[u], a.k, x, a.unit_kmers) for u, x in shared.items()}
        else:
            n_star = int(np.ceil(np.log(a.detect_miss) / np.log(1 - a.detect_aai**a.k)))
            unit, shared = detect_units(tree, n_star)
            chosen = {u: x[: 2 * n_star] for u, x in shared.items()}
        multi = np.bincount(unit, minlength=n)[unit] > 1
        own_k = np.arange(len(h)) < len(h) - n_ladder  # rows of the index's k, not the ladder's
        drop = (pos < n) & own_k & multi[np.minimum(pos, n - 1)]
        members = pl.DataFrame({"u": unit, "pos": np.arange(n)}).filter(pl.Series(multi))
        add = pl.DataFrame(
            {"u": np.repeat(list(chosen), [len(x) for x in chosen.values()]),
             "h": np.concatenate(list(chosen.values())) if chosen else np.zeros(0, U64)}
        ).join(members, on="u")  # fmt: skip
        pos = np.concatenate([pos[~drop], add["pos"].to_numpy()])
        h = np.concatenate([h[~drop], add["h"].to_numpy()])
    if a.unit_radius:  # every leaf of a unit holds the union of its members' k-mers
        unit = units(tree, a.unit_radius)
        held = pl.DataFrame({"u": unit[pos], "h": h}).unique()
        spread = held.join(pl.DataFrame({"u": unit, "pos": np.arange(n)}), on="u")
        pos, h = spread["pos"].to_numpy(), spread["h"].to_numpy()
    leaf = lay["order"][pos]
    o = np.lexsort((h, leaf))
    # every hash a leaf does not hold is known absent (no missing data above t_base)
    top = U64(ALL) if a.window or a.ladder else lo_base
    sets = p0.Sets(p0.csr(leaf[o], h[o], n), h[o], np.zeros(n, np.uint8), np.array([top]))
    ch = _core.tree_schemes(
        n, f.left, f.right, f.tops, sets.offsets.astype(U64), sets.hashes, sets.level,
        sets.level_hash, emit=True,
    )  # fmt: skip
    if a.cut_detect:  # leaves below a detectable node hold its state: drop the diffs below
        n_star = int(np.ceil(np.log(a.detect_miss) / np.log(1 - a.cut_detect**a.k)))
        keep_ch = ~below_cuts(f, ch["change_node"], ch["change_gain"], n_star)
        ch = {x: v[keep_ch] for x, v in ch.items() if x.startswith("change_")}
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
    return table, pieces, n + len(leafless), key_bits(table, pl.concat(info))


def spread_shared(seq: bytes, k: int, shared: NDArray[np.uint64], n_sel: int) -> NDArray[np.uint64]:
    """About ``n_sel`` of ``shared`` spread along ``seq``: the least shared hash in each
    window of n_k // (n_sel / 2) k-mers (robust winnowing restricted to shared k-mers)."""
    ks = [seq[i : i + k] for i in range(len(seq) - k + 1)]
    o = _core.hash_proteins(ks, k)
    none = np.iinfo(np.uint64).max
    h = np.full(len(ks), none, dtype=np.uint64)
    h[o["seq"]] = o["hash"]
    h[~np.isin(h, shared)] = none
    _, m = winnow(h, np.array([len(h)]), max(1, 2 * len(h) // max(n_sel, 1)))
    return np.unique(m[m != none])


def below_cuts(
    f: p0.Forest, node: NDArray[np.uint32], gain: NDArray[np.bool_], n_star: int
) -> NDArray[np.bool_]:
    """Per change, whether its node lies strictly below a cut. A cut is the highest node v
    whose every leaf shares at least ``n_star`` of v's Fitch state: |state(v)| less the
    losses on the path from v to the leaf (a k-mer gained and lost below v counted too, so
    the share is a lower bound). Leaves below a cut then hold its state."""
    n, m = f.n, len(f.left)
    total = n + m
    parent = np.full(total, -1, dtype=np.int64)
    parent[f.left] = np.arange(n, total)
    parent[f.right] = np.arange(n, total)
    node = node.astype(np.int64)
    g_at = np.bincount(node[gain], minlength=total)
    l_at = np.bincount(node[~gain], minlength=total)
    size = np.zeros(total, dtype=np.int64)
    lost = np.zeros(total, dtype=np.int64)
    # ponytail: Python loops over nodes (~1 s per 10^6); vectorised by level at scale
    for v in range(total - 1, -1, -1):  # parents have larger ids: top down
        p = parent[v]
        size[v] = (size[p] if p >= 0 else 0) + g_at[v] - l_at[v]
        lost[v] = (lost[p] if p >= 0 else 0) + l_at[v]
    most = lost.copy()  # most losses from the root to any leaf below
    for i in range(m):  # children have smaller ids: bottom up
        most[n + i] = max(most[f.left[i]], most[f.right[i]])
    ok = most - lost <= size - n_star
    cut = np.zeros(total, dtype=bool)
    stack = [int(t) for t in f.tops]
    while stack:
        v = stack.pop()
        if v >= n and ok[v]:
            cut[v] = True
        elif v >= n:
            stack += [int(f.left[v - n]), int(f.right[v - n])]
    inner = np.zeros(total, dtype=bool)
    for v in range(total - 1, -1, -1):
        p = parent[v]
        inner[v] = p >= 0 and (cut[p] or inner[p])
    return inner[node]


def soft_units(
    tree: survey.Tree, n_sel: int, tau: float, held: Callable[[int], NDArray[np.uint64]]
) -> tuple[NDArray[np.int64], dict[int, NDArray[np.uint64]]]:
    """Units on conserved k-mers: as :func:`detect_units`, but a unit's set is the k-mers held
    by at least ``tau`` of its members, not by all. A leaf joins its parent's unit while that
    set keeps at least ``n_sel`` k-mers and the leaf holds at least ``tau`` of it, so every
    member stays covered. Returns each position's unit and each multi-leaf unit's set."""
    f, n = tree.forest, tree.forest.n
    pos_of_leaf = np.empty(n, dtype=np.int64)
    pos_of_leaf[tree.lay["order"]] = np.arange(n)
    a = pos_of_leaf[np.concatenate([f.mst_a, f.mst_b])]
    b = pos_of_leaf[np.concatenate([f.mst_b, f.mst_a])]
    o = np.argsort(a, kind="stable")
    a, b = a[o], b[o]
    start = np.searchsorted(a, np.arange(n + 1))
    unit = np.full(n, -1, dtype=np.int64)
    counts: dict[int, tuple[NDArray[np.uint64], NDArray[np.int64], int]] = {}
    # ponytail: a Python walk merging count arrays; a Rust kernel at release scale
    for r in pos_of_leaf[f.root_nodes].tolist():
        unit[r] = r
        counts[r] = (held(r), np.ones(len(held(r)), np.int64), 1)
        stack = [r]
        while stack:
            v = stack.pop()
            for c in b[start[v] : start[v + 1]].tolist():
                if unit[c] >= 0:
                    continue
                u = unit[v]
                hs, cs, m = counts[u]
                hc = held(c)
                cur = hs[cs >= tau * m]
                allh = np.concatenate([hs, hc])
                keys, inv = np.unique(allh, return_inverse=True)
                new = np.bincount(inv, weights=np.concatenate([cs, np.ones(len(hc))])).astype(
                    np.int64
                )
                kept = keys[new >= tau * (m + 1)]
                if len(kept) >= n_sel and (np.isin(cur, hc).mean() if len(cur) else 1.0) >= tau:
                    unit[c], counts[u] = u, (keys, new, m + 1)
                else:
                    unit[c], counts[c] = c, (hc, np.ones(len(hc), np.int64), 1)
                stack.append(c)
    sizes = np.bincount(unit, minlength=n)
    return unit, {u: hs[cs >= tau * m] for u, (hs, cs, m) in counts.items() if sizes[u] > 1}


def detect_units(
    tree: survey.Tree, n_star: int, held: Callable[[int], NDArray[np.uint64]] | None = None
) -> tuple[NDArray[np.int64], dict[int, NDArray[np.uint64]]]:
    """Units by detectability: each spanning-tree component walked from its root, a leaf
    joining its parent's unit while the unit's members still share at least ``n_star``
    k-mers (hashes <= t_cap), enough to detect each member at the target AAI. Returns each
    position's unit (a position) and, per multi-leaf unit, its shared hashes. ``held``
    gives a position's set (default: its hashes <= t_cap)."""
    f, n = tree.forest, tree.forest.n
    pos_of_leaf = np.empty(n, dtype=np.int64)
    pos_of_leaf[tree.lay["order"]] = np.arange(n)
    a = pos_of_leaf[np.concatenate([f.mst_a, f.mst_b])]
    b = pos_of_leaf[np.concatenate([f.mst_b, f.mst_a])]
    o = np.argsort(a, kind="stable")
    a, b = a[o], b[o]
    start = np.searchsorted(a, np.arange(n + 1))
    if held is None:
        held = lambda v: tree.h[tree.row[v] : tree.row[v + 1]]  # noqa: E731
    unit = np.full(n, -1, dtype=np.int64)
    shared: dict[int, NDArray[np.uint64]] = {}
    # ponytail: a Python walk with set intersections; a Rust kernel at release scale
    for r in pos_of_leaf[f.root_nodes].tolist():
        unit[r] = r
        shared[r] = held(r)
        stack = [r]
        while stack:
            v = stack.pop()
            for c in b[start[v] : start[v + 1]].tolist():
                if unit[c] >= 0:
                    continue
                u = unit[v]
                both = np.intersect1d(shared[u], held(c), assume_unique=True)
                if len(both) >= n_star:
                    unit[c], shared[u] = u, both
                else:
                    unit[c], shared[c] = c, held(c)
                stack.append(c)
    sizes = np.bincount(unit, minlength=n)
    return unit, {u: x for u, x in shared.items() if sizes[u] > 1}


def units(tree: survey.Tree, radius: float) -> NDArray[np.int64]:
    """Each position's unit: spanning-tree components walked from their roots, a leaf joining
    its parent's unit while its summed identity distance from the unit's first leaf is at
    most ``radius``. An edge's distance is 1 - p, from its weight w = |A Δ B| of dense k-mer
    sets: p^k = 1 - w / (|A| + |B|), set sizes ~ length - k + 1."""
    f, n = tree.forest, tree.forest.n
    k = tree.args.topology_k
    size = np.empty(n, dtype=np.float64)
    size[tree.lay["order"]] = [max(len(tree.seqs[i]) - k + 1, 1) for i in tree.rep_by_pos]
    share = np.clip(1 - f.mst_w / (size[f.mst_a] + size[f.mst_b]), 0, 1)
    dist = 1 - share ** (1 / k)
    a = np.concatenate([f.mst_a, f.mst_b])
    b = np.concatenate([f.mst_b, f.mst_a])
    d = np.concatenate([dist, dist])
    o = np.argsort(a, kind="stable")
    a, b, d = a[o], b[o], d[o]
    start = np.searchsorted(a, np.arange(n + 1))
    unit = np.full(n, -1, dtype=np.int64)
    from_rep = np.zeros(n)
    # ponytail: a Python walk over leaves (~1 s per 10^6); a Rust kernel at release scale
    for r in f.root_nodes.tolist():
        unit[r] = r
        stack = [r]
        while stack:
            v = stack.pop()
            for c, dc in zip(
                b[start[v] : start[v + 1]].tolist(),
                d[start[v] : start[v + 1]].tolist(),
                strict=True,
            ):
                if unit[c] >= 0:
                    continue
                if from_rep[v] + dc <= radius:
                    unit[c], from_rep[c] = unit[v], from_rep[v] + dc
                else:
                    unit[c], from_rep[c] = c, 0.0
                stack.append(c)
    pos_of_leaf = np.empty(n, dtype=np.int64)
    pos_of_leaf[tree.lay["order"]] = np.arange(n)
    return pos_of_leaf[unit][tree.lay["order"]]


def key_bits(table: pl.DataFrame, info: pl.DataFrame) -> pl.DataFrame:
    """(hash, k, bits) per key: a hit's evidence, -log2 q0 - log2 g, g the key's gains (its
    positive postings: each an independent region of the tree a hit seeds)."""
    g = table.filter(pl.col("sign") > 0).group_by("hash").agg(g=pl.len())
    return (
        g.join(info.unique("hash"), on="hash", how="left")
        .with_columns(bits=pl.col("bits0").fill_null(0.0) - pl.col("g").cast(pl.Float64).log(2))
        .select("hash", pl.col("k").fill_null(0), "bits")
        .sort("hash")
    )


def residue_log_freq(seqs: list[bytes]) -> NDArray[np.float64]:
    """log2 of each byte's frequency over the sequences (a floor for unseen bytes)."""
    counts = np.bincount(np.frombuffer(b"".join(seqs), dtype=np.uint8), minlength=256) + 1e-9
    return np.log2(counts / counts.sum())


def kmer_info(seqs: list[bytes], k: int, log_f: NDArray[np.float64]) -> pl.DataFrame:
    """(seq, hash, bits0) for every valid k-mer: bits0 = -log2 q0, q0 its chance of occurring
    at a protein window from residue composition (low complexity: fewer bits)."""
    parts = []
    for start in range(0, len(seqs), 20_000):
        batch = seqs[start : start + 20_000]
        n_k = np.array([max(len(s) - k + 1, 0) for s in batch])
        # ponytail: each k-mer hashed as its own sequence; a positional kernel at scale
        o = _core.hash_proteins([s[i : i + k] for s in batch for i in range(len(s) - k + 1)], k)
        lf = np.cumsum(np.r_[0, log_f[np.frombuffer(b"".join(batch), dtype=np.uint8)]])
        first_res = np.cumsum([0] + [len(s) for s in batch])[:-1]
        within = np.arange(n_k.sum()) - np.repeat(np.cumsum(n_k) - n_k, n_k)
        at = np.repeat(first_res, n_k) + within  # each k-mer's first residue
        bits0 = -(lf[at + k] - lf[at])
        seq = np.repeat(np.arange(len(batch)), n_k) + start
        parts.append(pl.DataFrame({"seq": seq[o["seq"]], "hash": o["hash"],
                                   "bits0": bits0[o["seq"]]}))  # fmt: skip
    return pl.concat(parts).unique(["seq", "hash"])


def informative(info: pl.DataFrame, n_seq: int, rare: float) -> NDArray[np.int64]:
    """Distinct k-mers per sequence with bits0 >= ``rare``."""
    got = info.filter(pl.col("bits0") >= rare).group_by("seq").agg(pl.len())
    out = np.zeros(n_seq, dtype=np.int64)
    out[got["seq"].to_numpy()] = got["len"].to_numpy()
    return out


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
        help="k:window,...: a leaf with fewer than --min-kmers informative k-mers at the k "
        "before also holds window minimizers of k",
    )  # fmt: skip
    parser.add_argument("--min-kmers", type=int, default=64)
    parser.add_argument(
        "--detect-aai", type=float, default=0.0,
        help="units by detectability: members share enough k-mers to detect each at this AAI",
    )  # fmt: skip
    parser.add_argument("--detect-miss", type=float, default=0.007)
    parser.add_argument(
        "--cut-detect", type=float, default=0.0,
        help="drop the diffs below each highest node whose state still detects every leaf below "
        "at this AAI (n* shared k-mers, from --detect-miss); 0: off",
    )  # fmt: skip
    parser.add_argument(
        "--unit-tau",
        type=float,
        default=1.0,
        help="with --unit-kmers: a unit's set is the k-mers held by this share of its members",
    )
    parser.add_argument(
        "--unit-kmers", type=int, default=0,
        help="units on dense intersections: members share at least this many k-mers, and the "
        "unit indexes about this many of them, spread (0: off)",
    )  # fmt: skip
    parser.add_argument(
        "--unit-radius", type=float, default=0.0,
        help="leaves within this identity distance of their unit's first leaf share the "
        "union of their k-mers (0: a set per leaf)",
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
    table, pieces, n_pos, bits = postings(tree, a)
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
    np.save(out / "keys.npy", bits["hash"].to_numpy())
    np.save(out / "key_bits.npy", bits["bits"].to_numpy().astype(np.float32))
    np.save(out / "key_k.npy", bits["k"].to_numpy().astype(np.uint8))
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
