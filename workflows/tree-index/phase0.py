"""Tree index, phase 0: measure the k-mer tree on a sample of MGnify, no index
(docs/tree-index-plan.md, Implementation plan, phase 0).

One subcommand per pipeline step (``main.nf``):

- ``fixture``: synthetic protein families as members buckets (the test profile).
- ``sample``: one members bucket -> the sample's clusters (``cluster_rep % --sample == 0``)
  and their stored k-mers (``core_hashes.npy``), for ``stream``'s cross-cluster scan.
- ``concat``: the buckets' samples -> ``sample.parquet`` and ``core_hashes.npy``; or the
  buckets' representatives -> ``extra.parquet``.
- ``stream``: one bucket of the whole release -> counts by hash range (set digests,
  pigeonhole signatures, k-mer holders in hash slices) and the cluster representatives
  outside the sample that share at least ``--min-shared`` stored k-mers with it.
- ``reduce``: one hash range of every bucket -> group-size histograms.
- ``combine``: every range -> ``stream.tsv`` and ``stream.json``.
- ``measure``: the sample's trees -> ``measure.json`` and TSVs (see ``measure``).
- ``report``: ``measure.json`` and today's tier-1 ``meta.json`` -> ``gate.tsv``.

Stored sets follow the plan: S\\_v = k-mers of v with hash <= t\\_v, t\\_v = max(t\\_base,
n\\_min / (L - k + 1)) capped at t\\_cap and rounded up to one of ``--levels`` log-spaced
rates (Variable thresholds); ``--no-floors`` stores every sequence at t\\_base.
"""

import argparse
import json
import time
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

import numpy as np
import polars as pl
from numpy.typing import NDArray
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from kmer_functional_profiler import _core
from kmer_functional_profiler.index import _dense_ids, _gather, _mix, _set_hashes, mask_adapters

ROOT_SEED: Final = 0x9E3779B97F4A7C15  # part and signature hashes
BATCH: Final = 500_000  # sequences hashed at once
U64: Final = np.uint64


# ---------------------------------------------------------------------------- stored sets


@dataclass(frozen=True)
class Rates:
    k: int
    t_base: float
    n_min: int
    t_cap: float
    levels: int
    floors: bool = True

    @property
    def rates(self) -> NDArray[np.float64]:
        if self.levels == 1:
            return np.array([self.t_base])
        return self.t_base * (self.t_cap / self.t_base) ** (
            np.arange(self.levels) / (self.levels - 1)
        )

    @property
    def level_hash(self) -> NDArray[np.uint64]:
        return np.array([_core.max_hash(float(r)) for r in self.rates], dtype=U64)

    def level(self, lengths: NDArray[np.int64]) -> NDArray[np.uint8]:
        """Each sequence's level: the smallest rate at or above its t\\_v."""
        if not self.floors:
            return np.zeros(len(lengths), dtype=np.uint8)
        windows = np.maximum(lengths - self.k + 1, 1)
        t = np.clip(np.maximum(self.t_base, self.n_min / windows), None, self.t_cap)
        return np.searchsorted(self.rates, t * (1 - 1e-12)).astype(np.uint8)


def rates_of(args: argparse.Namespace, floors: bool = True) -> Rates:
    return Rates(args.k, args.t_base, args.n_min, args.t_cap, args.levels, floors)


@dataclass
class Sets:
    """CSR sets (sorted, distinct) with a level each; ``level_hash`` maps levels to hashes."""

    offsets: NDArray[np.int64]
    hashes: NDArray[np.uint64]
    level: NDArray[np.uint8]
    level_hash: NDArray[np.uint64]

    @property
    def n(self) -> int:
        return len(self.offsets) - 1

    @property
    def size(self) -> NDArray[np.int64]:
        return np.diff(self.offsets)

    @property
    def max_hash(self) -> NDArray[np.uint64]:
        return self.level_hash[self.level]

    def take(self, rows: NDArray[np.int64], level: NDArray[np.uint8] | None = None) -> "Sets":
        offsets, hashes = _gather(self.offsets, self.hashes, rows)
        return Sets(offsets, hashes, self.level[rows] if level is None else level, self.level_hash)

    def owner(self) -> NDArray[np.int64]:
        return np.repeat(np.arange(self.n), self.size)

    def kernel(self) -> tuple[NDArray[np.uint64], NDArray[np.uint64], NDArray[np.uint64]]:
        return self.offsets.astype(U64), self.hashes, self.max_hash


def hashed(
    seqs: list[bytes], k: int, max_hash: int, alphabet: str = "protein", min_hash: int = 0
) -> tuple[NDArray[np.int64], NDArray[np.uint64]]:
    """Distinct (sequence, hash) pairs with hash in ``min_hash..=max_hash``, sorted."""
    seq_parts, hash_parts = [], []
    for start in range(0, len(seqs), BATCH):
        out = _core.hash_proteins(
            seqs[start : start + BATCH], k, alphabet=alphabet, max_hash=max_hash, min_hash=min_hash
        )
        seq_parts.append(out["seq"].astype(np.int64) + start)
        hash_parts.append(out["hash"])
    seq = np.concatenate(seq_parts) if seq_parts else np.zeros(0, np.int64)
    h = np.concatenate(hash_parts) if hash_parts else np.zeros(0, U64)
    order = np.lexsort((h, seq))
    seq, h = seq[order], h[order]
    keep = np.ones(len(h), dtype=bool)
    keep[1:] = (seq[1:] != seq[:-1]) | (h[1:] != h[:-1])
    return seq[keep], h[keep]


def csr(seq: NDArray[np.int64], h: NDArray[np.uint64], n: int) -> NDArray[np.int64]:
    return np.concatenate([[0], np.cumsum(np.bincount(seq, minlength=n))]).astype(np.int64)


def stored_sets(seqs: list[bytes], rates: Rates) -> Sets:
    """Each sequence's stored set at its own level."""
    level = rates.level(np.array([len(s) for s in seqs], dtype=np.int64))
    level_hash = rates.level_hash
    seq, h = hashed(seqs, rates.k, int(level_hash.max()))
    keep = h <= level_hash[level[seq]]
    seq, h = seq[keep], h[keep]
    return Sets(csr(seq, h, len(seqs)), h, level, level_hash)


def dense_sets(seqs: list[bytes], k: int, rate: float = 1.0, alphabet: str = "protein") -> Sets:
    """Every sequence at one uniform rate (one level)."""
    level_hash = np.array([_core.max_hash(rate)], dtype=U64)
    seq, h = hashed(seqs, k, int(level_hash[0]), alphabet)
    return Sets(csr(seq, h, len(seqs)), h, np.zeros(len(seqs), np.uint8), level_hash)


def load_members(path: str | Path | list[str]) -> pl.DataFrame:
    """A members table (``protein_id``, ``cluster_rep``, ``full_length``, ``sequence``),
    adapters masked as the index build does, sorted by cluster."""
    members = pl.read_parquet(
        path, columns=["protein_id", "cluster_rep", "full_length", "sequence"]
    )
    return members.with_columns(mask_adapters(pl.col("sequence"))).sort("cluster_rep", "protein_id")


def seq_bytes(members: pl.DataFrame) -> list[bytes]:
    return members["sequence"].cast(pl.Binary).to_list()


def digests(sets: Sets) -> NDArray[np.uint64]:
    """Content hash of each set (rows of two 64-bit words); empty sets get (0, 0)."""
    out = np.zeros((sets.n, 2), dtype=U64)
    nonempty = sets.size > 0
    if nonempty.any():
        rows = np.flatnonzero(nonempty)
        sub = sets.take(rows)
        out[rows] = _set_hashes(sub.offsets, sub.hashes)
    return out


# ---------------------------------------------------------------------------- collapses


@dataclass
class Collapsed:
    leaves: Sets
    leaf_of: NDArray[np.int64]  # per sequence, its leaf (-1: empty set)
    stats: dict[str, Any]


def collapse(sets: Sets, full_length: NDArray[np.bool_]) -> tuple[Collapsed, dict[str, Any]]:
    """Equal-set and containment collapses (The structure, Collapses).

    Returns the leaves and the candidate pairs among distinct sets (reused by the build).
    A distinct set is contained when every k-mer it could share with another set (hashes at
    most both thresholds) is in it, and the other holds more there; it joins the closest
    such container, and chains resolve to the outermost.
    """
    nonempty = np.flatnonzero(sets.size > 0)
    ids, first = _dense_ids(digests(sets)[nonempty])
    n_distinct = len(first)
    level = np.zeros(n_distinct, dtype=np.uint8)
    np.maximum.at(level, ids, sets.level[nonempty])
    distinct = sets.take(nonempty[first], level)
    t0 = time.perf_counter()
    pairs = _core.overlap_pairs(*distinct.kernel())
    pair_s = time.perf_counter() - t0
    a, b = pairs["a"].astype(np.int64), pairs["b"].astype(np.int64)
    shared, weight = pairs["shared"].astype(np.int64), pairs["weight"].astype(np.int64)
    size, mh = distinct.size, distinct.max_hash
    a_low = mh[a] <= mh[b]  # the lower threshold's set is whole on the common range
    ra = np.where(a_low, size[a], weight + 2 * shared - size[b])
    rb = weight + 2 * shared - ra
    order_a = (size[a] < size[b]) | ((size[a] == size[b]) & (a < b))
    a_in_b = (ra == shared) & (rb > shared) & (shared > 0) & order_a
    b_in_a = (rb == shared) & (ra > shared) & (shared > 0) & ~order_a
    child = np.concatenate([a[a_in_b], b[b_in_a]])
    parent = np.concatenate([b[a_in_b], a[b_in_a]])
    w = np.concatenate([weight[a_in_b], weight[b_in_a]])
    container = np.arange(n_distinct)
    if len(child):
        best = (
            pl.DataFrame({"child": child, "parent": parent, "w": w})
            .sort("child", "w", "parent")
            .unique("child", keep="first")
        )
        container[best["child"].to_numpy()] = best["parent"].to_numpy()
        for _ in range(64):  # pointer jumping: containers grow, so no cycles
            nxt = container[container]
            if (nxt == container).all():
                break
            container = nxt
    is_leaf = container == np.arange(n_distinct)
    leaf_index = np.cumsum(is_leaf) - 1
    leaves = distinct.take(np.flatnonzero(is_leaf))
    leaf_of = np.full(sets.n, -1, dtype=np.int64)
    leaf_of[nonempty] = leaf_index[container[ids]]
    contained = np.zeros(sets.n, dtype=bool)
    contained[nonempty] = ~is_leaf[ids]
    partial = ~full_length
    stats = {
        "sequences": sets.n,
        "empty_sets": int(sets.n - len(nonempty)),
        "distinct_sets": n_distinct,
        "equal_set_collapsed": int(len(nonempty) - n_distinct),
        "contained_sets": int((~is_leaf).sum()),
        "leaves": leaves.n,
        "contained_share_partial": float(contained[partial].mean()) if partial.any() else None,
        "contained_share_full": float(contained[full_length].mean()) if full_length.any() else None,
        "stored_kmers": int(sets.size.sum()),
        "pair_join_s": pair_s,
        "pair_probes": int(pairs["probes"]),
    }
    keep = is_leaf[a] & is_leaf[b]
    leaf_pairs = {
        "a": leaf_index[a[keep]],
        "b": leaf_index[b[keep]],
        "weight": weight[keep],
        "shared": shared[keep],
    }
    return Collapsed(leaves, leaf_of, stats), leaf_pairs


# ---------------------------------------------------------------------------- trees


@dataclass
class Forest:
    n: int
    left: NDArray[np.uint32]
    right: NDArray[np.uint32]
    level: NDArray[np.uint32]
    tops: NDArray[np.uint32]
    mst_a: NDArray[np.int64]  # internal spanning-tree edges (leaf, leaf, weight)
    mst_b: NDArray[np.int64]
    mst_w: NDArray[np.int64]
    weight: int  # total, root edges included
    root_nodes: NDArray[np.int64]  # nodes whose spanning-tree edge goes to the empty root


def spanning_forest(
    n: int,
    a: NDArray[np.int64],
    b: NDArray[np.int64],
    w: NDArray[np.int64],
    root: NDArray[np.int64] | None,
) -> Forest:
    """Kruskal over (a, b, w), plus an edge from each node to the virtual empty root of
    weight ``root`` (its stored size) when given; then the merge forest of the tree's edges
    between nodes (the root's edges become the forest's cuts)."""
    u32 = np.uint32
    if root is not None:
        a = np.concatenate([a, np.arange(n)])
        b = np.concatenate([b, np.full(n, n)])
        w = np.concatenate([w, root])
    full = _core.kruskal(n + (root is not None), a.astype(u32), b.astype(u32), w.astype(u32))
    mst = full["mst"]
    inner = mst & (b < n) if root is not None else mst
    forest = _core.kruskal(n, a[inner].astype(u32), b[inner].astype(u32), w[inner].astype(u32))
    return Forest(
        n,
        forest["left"],
        forest["right"],
        forest["level"],
        forest["tops"],
        a[inner],
        b[inner],
        w[inner],
        int(w[mst].sum()),
        a[mst & (b == n)] if root is not None else np.zeros(0, np.int64),
    )


def schemes(forest: Forest, sets: Sets, tops: NDArray[np.uint32] | None = None) -> pl.DataFrame:
    """Per k-mer postings of each diff scheme on the forest (``tops`` overrides its roots)."""
    out = _core.tree_schemes(
        forest.n,
        forest.left,
        forest.right,
        forest.tops if tops is None else tops,
        sets.offsets.astype(U64),
        sets.hashes,
        sets.level,
        sets.level_hash,
    )
    return pl.DataFrame(out)


def node_sizes(
    forest: Forest, sets: Sets
) -> tuple[NDArray[np.int64], dict[str, NDArray[np.uint32]]]:
    """Stored k-mers under each node, and the forest's layout."""
    lay = _core.tree_layout(forest.n, forest.left, forest.right, forest.tops)
    prefix = np.concatenate([[0], np.cumsum(sets.size[lay["order"]])])
    lo, hi = lay["lo"].astype(np.int64), lay["hi"].astype(np.int64)
    reached = lay["lo"] != np.iinfo(np.uint32).max
    size = np.zeros(len(lo), dtype=np.int64)
    size[reached] = prefix[hi[reached]] - prefix[lo[reached]]
    return size, lay


def chunk_roots(forest: Forest, size: NDArray[np.int64], budget: float) -> NDArray[np.uint32]:
    """Top-down: a subtree of stored size <= ``budget`` (or a leaf) is a chunk; a larger one
    splits into its children (Chunks by similarity)."""
    roots, stack = [], [int(t) for t in forest.tops]
    while stack:
        x = stack.pop()
        if x < forest.n or size[x] <= budget:
            roots.append(x)
        else:
            i = x - forest.n
            stack += [int(forest.left[i]), int(forest.right[i])]
    return np.array(sorted(roots), dtype=np.uint32)


# ---------------------------------------------------------------------------- staged build


def signatures(sets: Sets, parts: int) -> pl.DataFrame:
    """Pigeonhole signatures: per set and part (a second hash of the k-mer modulo
    ``parts``), a digest of the set's k-mers in that part; empty parts emit nothing."""
    owner = sets.owner()
    part = (_mix(sets.hashes ^ U64(ROOT_SEED)) % U64(parts)).astype(np.int64)
    order = np.lexsort((sets.hashes, part, owner))
    owner, part, h = owner[order], part[order], sets.hashes[order]
    key = owner * parts + part
    starts = np.flatnonzero(np.concatenate([[True], key[1:] != key[:-1]]))
    offsets = np.concatenate([starts, [len(key)]]).astype(np.int64)
    dig = _set_hashes(offsets, h)
    return pl.DataFrame(
        {"node": owner[starts], "part": part[starts], "d0": dig[:, 0], "d1": dig[:, 1]}
    )


def group_pairs(
    groups: pl.DataFrame, cap: int, order_col: str
) -> tuple[pl.DataFrame, dict[str, int]]:
    """Candidate pairs within groups (``group``, ``node``, ``order_col``, ``size``): all
    pairs in a group of at most ``cap`` nodes; above it, each node is linked to its sort
    neighbours by ``order_col`` and to the group's largest node (Linclust's centre)."""
    sized = groups.with_columns(g=pl.len().over("group"))
    small = sized.filter(pl.col("g") <= cap).select("group", "node")
    pairs = [
        small.join(small, on="group", suffix="_b")
        .filter(pl.col("node") < pl.col("node_b"))
        .select(a="node", b="node_b")
    ]
    large = sized.filter(pl.col("g") > cap).sort("group", order_col, "node")
    if large.height:
        nxt = large.select(a="node", b=pl.col("node").shift(-1).over("group")).drop_nulls()
        centre = (
            large.sort("group", "size", "node", descending=[False, True, False])
            .group_by("group", maintain_order=True)
            .first()
            .select("group", c="node")
        )
        star = large.join(centre, on="group").select(a="node", b="c")
        pairs += [nxt, star]
    out = (
        pl.concat(pairs)
        .filter(pl.col("a") != pl.col("b"))
        .select(a=pl.min_horizontal("a", "b"), b=pl.max_horizontal("a", "b"))
        .unique()
    )
    sizes = sized.group_by("group").agg(g=pl.first("g"))["g"]
    stats = {
        "groups": sizes.len(),
        "largest_group": int(sizes.max() or 0),  # type: ignore[arg-type]
        "capped_groups": int((sizes > cap).sum()),
        "candidates": out.height,
    }
    return out, stats


def weigh(sets: Sets, pairs: pl.DataFrame) -> pl.DataFrame:
    a = pairs["a"].to_numpy().astype(np.uint32)
    b = pairs["b"].to_numpy().astype(np.uint32)
    got = _core.pair_weights(*sets.kernel(), a, b)
    return pairs.with_columns(shared=got["shared"], weight=got["weight"])


def staged_build(
    leaves: Sets, args: argparse.Namespace
) -> tuple[Forest, dict[str, Any], pl.DataFrame]:
    """Stages 3-4 of Building the tree (stage 2, equal sets, is the collapse): pigeonhole
    signatures for weights <= d, then a rarity-ordered prefix join on one representative per
    component, with sorted-neighbourhood candidates for capped postings; Kruskal over both
    plus the root's edges. Returns the forest, its stats and the prefix posting lengths."""
    d, size = args.d, leaves.size
    stats: dict[str, Any] = {}
    t0 = time.perf_counter()
    sig = signatures(leaves, d + 1)
    node = sig["node"].to_numpy()
    groups = sig.with_columns(
        group=pl.struct("part", "d0", "d1").hash(), size=size[node], full=digests(leaves)[node, 0]
    )
    cand3, s3 = group_pairs(groups.select("group", "node", "full", "size"), args.group_cap, "full")
    cand3 = weigh(leaves, cand3)
    edges3 = cand3.filter(pl.col("weight") <= d)
    stats |= {f"stage3_{k}": v for k, v in s3.items()} | {
        "stage3_false_positives": cand3.height - edges3.height,
        "stage3_edges": edges3.height,
    }
    n = leaves.n
    graph = coo_matrix(
        (np.ones(edges3.height), (edges3["a"].to_numpy(), edges3["b"].to_numpy())), shape=(n, n)
    )
    _, comp = connected_components(graph, directed=False)
    reps = (
        pl.DataFrame({"node": np.arange(n), "comp": comp, "size": size})
        .sort("comp", "size", "node", descending=[False, True, False])
        .group_by("comp", maintain_order=True)
        .first()["node"]
        .to_numpy()
    )
    stats["stage4_representatives"] = len(reps)
    # rarity order among representatives, rarest first; prefix of floor(|S|/2) + 1
    owner = np.repeat(reps, size[reps])
    rep_sets = leaves.take(reps)
    inv = pl.DataFrame({"node": owner, "hash": rep_sets.hashes, "size": size[owner]})
    inv = inv.with_columns(holders=pl.len().over("hash"))
    prefix = (
        inv.sort("node", "holders", "hash")
        .with_columns(rank=pl.int_range(pl.len()).over("node"))
        .filter(pl.col("rank") < pl.col("size") // 2 + 1)
    )
    capped = prefix.filter(pl.col("holders") > args.posting_cap)
    probe = prefix.filter(pl.col("holders") <= args.posting_cap).select("node", "hash", "size")
    cand4 = (
        probe.join(inv.select("hash", node_b="node", size_b="size"), on="hash")
        .filter(
            (pl.col("size") < pl.col("size_b"))
            | ((pl.col("size") == pl.col("size_b")) & (pl.col("node") < pl.col("node_b")))
        )
        .select(a=pl.min_horizontal("node", "node_b"), b=pl.max_horizontal("node", "node_b"))
        .unique()
    )
    stats["stage4_prefix_candidates"] = cand4.height
    stats["stage4_capped_postings"] = capped["hash"].n_unique()
    # sorted neighbourhood for representatives with a capped prefix posting
    capped_reps = capped["node"].unique().to_numpy()
    sn = []
    if len(capped_reps):
        for r in range(args.sn_orders):
            key = _mix(rep_sets.hashes ^ U64(ROOT_SEED + r + 1))
            bottom = pl.DataFrame({"node": owner, "key": key}).group_by("node").agg(pl.min("key"))
            ranked = bottom.sort("key", "node").with_columns(pos=pl.int_range(pl.len()))
            mine = ranked.filter(pl.col("node").is_in(capped_reps.tolist()))
            for off in range(1, args.sn_window + 1):
                for sign in (1, -1):
                    nb = mine.select("node", pos=pl.col("pos") + sign * off).join(
                        ranked.select("pos", node_b="node"), on="pos"
                    )
                    sn.append(nb.select(a="node", b="node_b"))
        snp = (
            pl.concat(sn)
            .select(a=pl.min_horizontal("a", "b"), b=pl.max_horizontal("a", "b"))
            .unique()
        )
        cand4 = pl.concat([cand4, snp]).unique()
        stats["stage4_neighbourhood_candidates"] = snp.height
    cand4 = weigh(leaves, cand4).filter(pl.col("shared") > 0)
    larger = np.maximum(size[cand4["a"].to_numpy()], size[cand4["b"].to_numpy()])
    edges4 = cand4.filter(pl.col("weight") <= pl.Series(larger))
    stats["stage4_edges"] = edges4.height
    edges = pl.concat([edges3.select("a", "b", "weight"), edges4.select("a", "b", "weight")])
    forest = spanning_forest(
        n,
        edges["a"].to_numpy().astype(np.int64),
        edges["b"].to_numpy().astype(np.int64),
        edges["weight"].to_numpy().astype(np.int64),
        size,
    )
    stats["staged_build_s"] = time.perf_counter() - t0
    stats["staged_weight"] = forest.weight
    lengths = prefix.group_by("holders").len().rename({"holders": "posting", "len": "n"})
    return forest, stats, lengths.sort("posting")


# ---------------------------------------------------------------------------- measure


def histogram(values: NDArray[Any], name: str) -> pl.DataFrame:
    v, c = np.unique(values, return_counts=True)
    return pl.DataFrame({name: v, "n": c})


def sequence_mst_postings(forest: Forest, sets: Sets) -> int:
    """The sequence-node scheme's postings on the spanning tree: each node's state is its
    stored set, so an edge costs the whole symmetric difference (a parent's k-mers above the
    child's threshold are deleted too), and a root edge the node's set."""
    every = np.full(sets.n, np.iinfo(np.uint64).max, dtype=U64)
    got = _core.pair_weights(
        sets.offsets.astype(U64),
        sets.hashes,
        every,
        forest.mst_a.astype(np.uint32),
        forest.mst_b.astype(np.uint32),
    )
    return int(got["weight"].sum()) + int(sets.size[forest.root_nodes].sum())


def tree_postings(forest: Forest, sets: Sets) -> tuple[pl.DataFrame, dict[str, Any]]:
    t0 = time.perf_counter()
    per_kmer = schemes(forest, sets)
    fitch_s = time.perf_counter() - t0
    t0 = time.perf_counter()
    sequence_mst = sequence_mst_postings(forest, sets)
    sequence_mst_s = time.perf_counter() - t0
    totals = {
        "fitch": int(per_kmer["fitch"].sum()),
        "majority": int(per_kmer["majority"].sum()),
        "add_only": int(per_kmer["add_only"].sum()),
        "sequence_mst": sequence_mst,
        "spanning_tree_weight": forest.weight,  # edges over hashes both ends could hold
        "flat": int(per_kmer["holders"].sum()),
        "distinct_kmers": per_kmer.height,
        "fitch_s": fitch_s,
        "sequence_mst_s": sequence_mst_s,
        "tops": len(forest.tops),
        "internal_nodes": len(forest.left),
    }
    return per_kmer, totals


def measure_variant(
    members: pl.DataFrame,
    seqs: list[bytes],
    rates: Rates,
    out: Path,
    tag: str,
    args: argparse.Namespace,
) -> tuple[dict[str, Any], Collapsed]:
    """Collapses, the exact tree, the schemes and (for the main variant) the staged build."""
    sets = stored_sets(seqs, rates)
    col, pairs = collapse(sets, members["full_length"].to_numpy())
    leaves = col.leaves
    t0 = time.perf_counter()
    exact = spanning_forest(leaves.n, pairs["a"], pairs["b"], pairs["weight"], leaves.size)
    kruskal_s = time.perf_counter() - t0
    per_kmer, totals = tree_postings(exact, leaves)
    result: dict[str, Any] = {
        "collapse": col.stats,
        "exact_weight": exact.weight,
        "exact_candidates": len(pairs["a"]),
        "exact_kruskal_s": kruskal_s,
        "postings": totals,
    }
    singles = exact.tops[exact.tops < leaves.n].astype(np.int64)
    base = leaves.take(singles)
    at_base = int((base.hashes <= rates.level_hash[0]).sum())
    result["postings_floors_on_components_only"] = totals["fitch"] - (
        int(base.size.sum()) - at_base
    )
    result["singleton_leaves"] = len(singles)
    if tag != "floors":
        return result, col
    # staged build against the exact tree
    staged, stage_stats, prefix = staged_build(leaves, args)
    prefix.write_csv(out / "prefix_postings.tsv", separator="\t")
    result["staged"] = stage_stats | {
        "weight_over_exact": staged.weight / max(exact.weight, 1),
        "staged_fitch": int(schemes(staged, leaves)["fitch"].sum()),
    }
    # per-k-mer counts, CI and cuts
    units = pl.DataFrame(
        {"hash": sets.hashes, "unit": members["cluster_rep"].to_numpy()[sets.owner()]}
    )
    n_units = units.group_by("hash").agg(units=pl.col("unit").n_unique())
    per_kmer = per_kmer.join(n_units, on="hash", how="left")
    per_kmer.write_parquet(out / "kmers.parquet")
    histogram(per_kmer["fitch"].to_numpy(), "changes").write_csv(
        out / "changes.tsv", separator="\t"
    )
    histogram(per_kmer["gains"].to_numpy(), "gains").write_csv(out / "gains.tsv", separator="\t")
    cuts = [
        {
            "cut": f"gains>{g}",
            "kmers_dropped": int((per_kmer["gains"] > g).sum()),
            "postings_saved": int(per_kmer.filter(pl.col("gains") > g)["fitch"].sum()),
        }
        for g in (2, 4, 8)
    ] + [
        {
            "cut": f"units>{args.max_groups}",
            "kmers_dropped": int((per_kmer["units"] > args.max_groups).sum()),
            "postings_saved": int(
                per_kmer.filter(pl.col("units") > args.max_groups)["fitch"].sum()
            ),
        }
    ]
    pl.DataFrame(cuts).write_csv(out / "cuts.tsv", separator="\t")
    result["cuts"] = cuts
    # cross-cluster sharing (a leaf's cluster: its first member's)
    leaf_unit = np.full(leaves.n, -1, dtype=np.int64)
    has = col.leaf_of >= 0
    reps = members["cluster_rep"].to_numpy()
    leaf_unit[col.leaf_of[has][::-1]] = reps[has][::-1]
    cross = leaf_unit[exact.mst_a] != leaf_unit[exact.mst_b]
    result["mst_edges"] = len(exact.mst_a)
    result["mst_edges_across_units"] = int(cross.sum())
    result["units"] = int(members["cluster_rep"].n_unique())
    # edge cost by identity
    identity_table(members, seqs, col, leaves, exact, pairs, args).write_csv(
        out / "identity.tsv", separator="\t"
    )
    # chunking at the stored rate
    result["chunks_stored"] = chunking(exact, leaves, totals["fitch"], args)
    result["density_stored"] = {"distinct": totals["distinct_kmers"], "total": totals["flat"]}
    np.save(out / "leaf_of.npy", col.leaf_of)
    save_forest(out / "forest.npz", exact)
    return result, col


def save_forest(path: Path, forest: Forest) -> None:
    np.savez(
        path,
        n=forest.n,
        left=forest.left,
        right=forest.right,
        level=forest.level,
        tops=forest.tops,
        mst_a=forest.mst_a,
        mst_b=forest.mst_b,
        mst_w=forest.mst_w,
        weight=forest.weight,
        root_nodes=forest.root_nodes,
    )


def chunking(
    forest: Forest, sets: Sets, base: int, args: argparse.Namespace
) -> list[dict[str, Any]]:
    """Duplicated postings and chunk counts at each budget (stored k-mers per chunk)."""
    size, _ = node_sizes(forest, sets)
    rows = []
    for budget in args.chunk_sizes:
        roots = chunk_roots(forest, size, budget)
        dup = int(schemes(forest, sets, roots)["fitch"].sum()) - base
        rows.append(
            {
                "budget": budget,
                "chunk_roots": len(roots),
                "duplicated_postings": dup,
                "duplicated_share": dup / max(base, 1),
                "components_over_budget": float((size[forest.tops] > budget).mean()),
            }
        )
    return rows


def identity_table(
    members: pl.DataFrame,
    seqs: list[bytes],
    col: Collapsed,
    leaves: Sets,
    forest: Forest,
    pairs: dict[str, NDArray[np.int64]],
    args: argparse.Namespace,
) -> pl.DataFrame:
    """Edge cost (weight over the larger stored set, so 1 = as costly as storing fresh)
    against overlap-alignment identity, for spanning-tree edges and for pairs sharing any
    stored k-mer (Recall-first: does the ~94% cut appear?)."""
    rng = np.random.default_rng(args.seed)
    lengths = np.array([len(s) for s in seqs])
    has = np.flatnonzero(col.leaf_of >= 0)
    # each leaf's longest member stands for it
    rep = (
        pl.DataFrame({"seq": has, "leaf": col.leaf_of[has], "len": lengths[has]})
        .sort("leaf", "len", "seq", descending=[False, True, False])
        .group_by("leaf", maintain_order=True)
        .first()
        .sort("leaf")["seq"]
        .to_numpy()
    )
    # pairs sharing any stored k-mer, among the leaves of a subsample of whole clusters
    near = np.unique(col.leaf_of[subsample(members, args.dense_seqs)])
    near = near[near >= 0]
    relaxed = _core.overlap_pairs(
        *leaves.take(near).kernel(), storage=False, max_holders=args.posting_cap
    )
    kinds = {
        "mst": (forest.mst_a, forest.mst_b, forest.mst_w),
        "shared_kmer": (
            near[relaxed["a"].astype(np.int64)],
            near[relaxed["b"].astype(np.int64)],
            relaxed["weight"].astype(np.int64),
        ),
    }
    frames = []
    size = leaves.size
    for kind, (a, b, w) in kinds.items():
        if len(a) > args.identity_pairs:
            pick = rng.choice(len(a), args.identity_pairs, replace=False)
            a, b, w = a[pick], b[pick], w[pick]
        got = _core.align_identity([seqs[i] for i in rep[a]], [seqs[i] for i in rep[b]])
        frames.append(
            pl.DataFrame(
                {
                    "kind": kind,
                    "identity": got["identity"],
                    "cost": w / np.maximum(np.maximum(size[a], size[b]), 1),
                }
            )
        )
    return (
        pl.concat(frames)
        .with_columns(bin=((pl.col("identity").cast(pl.Float64) * 100).floor() / 100).round(2))
        .group_by("kind", "bin")
        .agg(
            n=pl.len(),
            cost_mean=pl.col("cost").mean(),
            cost_median=pl.col("cost").median(),
            above_fresh=(pl.col("cost") > 1).mean(),
        )
        .sort("kind", "bin")
    )


def density(seqs: list[bytes], sets: Sets, k: int) -> list[dict[str, Any]]:
    """Distinct over total k-mers at uniform rates, and at the stored rates."""
    rows = []
    for rate in (0.06, 0.2, 1.0):
        distinct = _core.DistinctHashes(_core.max_hash(rate))
        total = 0
        for start in range(0, len(seqs), BATCH):
            seq, h = hashed(seqs[start : start + BATCH], k, _core.max_hash(rate))
            distinct.add(h)
            total += len(h)
        n_distinct = sum(len(distinct.take(i)) for i in range(len(distinct)))
        rows.append({"rate": rate, "distinct": n_distinct, "total": total})
    rows.append(
        {"rate": "stored", "distinct": len(np.unique(sets.hashes)), "total": len(sets.hashes)}
    )
    for r in rows:
        r["distinct_over_total"] = r["distinct"] / max(r["total"], 1)
    return rows


def subsample(members: pl.DataFrame, n_seqs: int) -> NDArray[np.int64]:
    """Rows of the first whole clusters (in ``cluster_rep`` order) up to ``n_seqs``."""
    sizes = members.group_by("cluster_rep", maintain_order=True).len()
    keep = sizes.filter(pl.col("len").cum_sum() <= n_seqs)["cluster_rep"]
    return np.flatnonzero(members["cluster_rep"].is_in(keep.implode()).to_numpy())


def one_topology(
    members: pl.DataFrame,
    seqs: list[bytes],
    stored: Path,
    rows: NDArray[np.int64],
    args: argparse.Namespace,
) -> list[dict[str, Any]]:
    """Diffs at t 1 (k 11) and at k 8 along the stored-rate tree, against a tree built at
    those settings, on a subsample (Density, one topology many densities).

    The stored tree's edges are reused with the dense weights between their leaves' longest
    members, plus dense pairs inside each stored leaf; Kruskal over those grafts every
    sequence onto the stored topology."""
    leaf_of = np.load(stored / "leaf_of.npy")[rows]
    tree = np.load(stored / "forest.npz")
    sub = [seqs[i] for i in rows]
    lengths = np.array([len(s) for s in sub])
    out = []
    for k in (11, 8):
        dense = dense_sets(sub, k)
        col, pairs = collapse(dense, members["full_length"].to_numpy()[rows])
        leaves = col.leaves
        own = spanning_forest(leaves.n, pairs["a"], pairs["b"], pairs["weight"], leaves.size)
        own_fitch = int(schemes(own, leaves)["fitch"].sum())
        # each dense leaf's stored leaf (its first member's), and each stored leaf's longest
        has = col.leaf_of >= 0
        dense_to_stored = np.full(leaves.n, -1, dtype=np.int64)
        dense_to_stored[col.leaf_of[has][::-1]] = leaf_of[has][::-1]
        longest = (
            pl.DataFrame({"stored": leaf_of[has], "dense": col.leaf_of[has], "len": lengths[has]})
            .sort("stored", "len", descending=[False, True])
            .group_by("stored", maintain_order=True)
            .first()
        )
        rep = dict(zip(longest["stored"].to_list(), longest["dense"].to_list(), strict=True))
        within = dense_to_stored[pairs["a"]] == dense_to_stored[pairs["b"]]
        between = [
            (rep[x], rep[y])
            for x, y in zip(tree["mst_a"].tolist(), tree["mst_b"].tolist(), strict=True)
            if x in rep and y in rep and rep[x] != rep[y]
        ]
        link_a = np.array([p[0] for p in between], dtype=np.int64)
        link_b = np.array([p[1] for p in between], dtype=np.int64)
        link_w = _core.pair_weights(
            *leaves.kernel(), link_a.astype(np.uint32), link_b.astype(np.uint32)
        )
        grafted = spanning_forest(
            leaves.n,
            np.concatenate([pairs["a"][within], link_a]),
            np.concatenate([pairs["b"][within], link_b]),
            np.concatenate([pairs["weight"][within], link_w["weight"].astype(np.int64)]),
            leaves.size,
        )
        out.append(
            {
                "k": k,
                "rate": 1.0,
                "sequences": len(rows),
                "leaves": leaves.n,
                "own_tree_fitch": own_fitch,
                "stored_topology_fitch": int(schemes(grafted, leaves)["fitch"].sum()),
                "own_tree_weight": own.weight,
                "stored_topology_weight": grafted.weight,
                "flat": int(leaves.size.sum()),
            }
        )
        if k == 11:
            out[-1]["chunks_dense"] = chunking(own, leaves, own_fitch, args)
    return out


def divergence(shared: NDArray[Any], weight: NDArray[Any], k: int = 6) -> NDArray[np.float64]:
    """1 - identity estimated from k-mer Jaccard J: shared share c = 2J / (1 + J) ~ p^k.
    Monotone in J, so the spanning trees are Jaccard's; in identity units for the radii."""
    union = weight.astype(np.float64) + shared
    jaccard = np.where(union > 0, shared / np.maximum(union, 1), 0.0)
    return 1 - (2 * jaccard / (1 + jaccard)) ** (1 / k)


def sensitive_topology(
    seqs: list[bytes],
    stored: Path,
    leaf_of_all: NDArray[np.int64],
    rows: NDArray[np.int64],
    leaves: Sets,
    args: argparse.Namespace,
) -> dict[str, Any]:
    """The exact k 6 Jaccard MST against the seeded build, on one representative per storage
    forest component of a subsample (A sensitive topology at MGnify scale)."""
    tree = np.load(stored / "forest.npz")
    forest = Forest(
        int(tree["n"]),
        tree["left"],
        tree["right"],
        tree["level"],
        tree["tops"],
        tree["mst_a"],
        tree["mst_b"],
        tree["mst_w"],
        int(tree["weight"]),
        tree["root_nodes"],
    )
    lay = _core.tree_layout(forest.n, forest.left, forest.right, forest.tops)
    comp = np.empty(forest.n, dtype=np.int64)
    for c, t in enumerate(forest.tops):
        comp[lay["order"][lay["lo"][t] : lay["hi"][t]]] = c
    lengths = np.array([len(s) for s in seqs])
    has = rows[leaf_of_all[rows] >= 0]
    picked = (
        pl.DataFrame({"seq": has, "comp": comp[leaf_of_all[has]], "len": lengths[has]})
        .sort("comp", "len", "seq", descending=[False, True, False])
        .group_by("comp", maintain_order=True)
        .first()
    )
    rep_seq = picked["seq"].to_numpy()
    rep_comp = picked["comp"].to_numpy()
    sub = [seqs[i] for i in rep_seq]
    n = len(sub)
    q = 1e6  # divergence in millionths, for Kruskal's integer weights
    dense = dense_sets(sub, 6)
    ex = _core.overlap_pairs(*dense.kernel(), storage=False)
    ea, eb = ex["a"].astype(np.int64), ex["b"].astype(np.int64)
    ed = divergence(ex["shared"], ex["weight"])
    t0 = time.perf_counter()
    exact = _core.kruskal(
        n, ea.astype(np.uint32), eb.astype(np.uint32), np.round(ed * q).astype(np.uint32)
    )
    exact_s = time.perf_counter() - t0
    # seeds: the m lowest reduced-alphabet k-mers; centre (longest) and sort neighbours
    t0 = time.perf_counter()
    s_seq, s_hash = hashed(sub, args.seed_k, _core.max_hash(1.0), alphabet=args.seed_alphabet)
    seeds = (
        pl.DataFrame({"node": s_seq, "hash": s_hash})
        .sort("node", "hash")
        .filter(pl.int_range(pl.len()).over("node") < args.seeds)
    )
    seeds = seeds.with_columns(size=np.array([len(s) for s in sub])[seeds["node"].to_numpy()])
    cand, sstats = group_pairs(
        seeds.select(group="hash", node="node", size="size", full="size"), 2, "full"
    )
    sketch = dense_sets(sub, 6, args.sketch_rate)
    sw = _core.pair_weights(
        *sketch.kernel(),
        cand["a"].to_numpy().astype(np.uint32),
        cand["b"].to_numpy().astype(np.uint32),
    )
    xw = _core.pair_weights(
        *dense.kernel(),
        cand["a"].to_numpy().astype(np.uint32),
        cand["b"].to_numpy().astype(np.uint32),
    )
    ca, cb = cand["a"].to_numpy().astype(np.int64), cand["b"].to_numpy().astype(np.int64)
    sd = divergence(sw["shared"], sw["weight"])
    cd = divergence(xw["shared"], xw["weight"])
    seeded = _core.kruskal(
        n, ca.astype(np.uint32), cb.astype(np.uint32), np.round(sd * q).astype(np.uint32)
    )
    seeded_s = time.perf_counter() - t0
    exact_edges = set(zip(ea[exact["mst"]].tolist(), eb[exact["mst"]].tolist(), strict=True))
    cand_edges = set(zip(ca.tolist(), cb.tolist(), strict=True))
    unreached = {}
    for radius in (0.05, 0.1, 0.15, 0.2, 0.3):
        near = ed <= radius
        for weights, name in ((cd, "exact"), (sd, "sketch")):
            keep = weights <= radius
            g = coo_matrix((np.ones(keep.sum()), (ca[keep], cb[keep])), shape=(n, n))
            _, lab = connected_components(g, directed=False)
            miss = lab[ea[near]] != lab[eb[near]]
            unreached[f"R{radius}_{name}_weights"] = float(miss.mean()) if near.any() else None
        unreached[f"R{radius}_pairs"] = int(near.sum())
    # storage (k 11 Fitch) along the sensitive topology, against the storage tree
    comp_leaf = (
        pl.DataFrame(
            {"comp": comp[leaf_of_all[has]], "leaf": leaf_of_all[has], "len": lengths[has]}
        )
        .sort("comp", "len", descending=[False, True])
        .group_by("comp", maintain_order=True)
        .first()
    )
    comp_to_leaf = dict(zip(comp_leaf["comp"].to_list(), comp_leaf["leaf"].to_list(), strict=True))
    in_sub = np.zeros(forest.n, dtype=bool)
    in_sub[np.unique(leaf_of_all[has])] = True
    keep = in_sub[forest.mst_a] & in_sub[forest.mst_b]
    sa = np.array([comp_to_leaf[rep_comp[i]] for i in ca[seeded["mst"]]], dtype=np.int64)
    sb = np.array([comp_to_leaf[rep_comp[i]] for i in cb[seeded["mst"]]], dtype=np.int64)
    offset = int(forest.mst_w.max(initial=0)) + 1  # storage edges merge first
    sens = spanning_forest(
        forest.n,
        np.concatenate([forest.mst_a[keep], sa]),
        np.concatenate([forest.mst_b[keep], sb]),
        np.concatenate(
            [forest.mst_w[keep], offset + np.round(sd[seeded["mst"]] * 1000).astype(np.int64)]
        ),
        None,
    )
    storage = spanning_forest(
        forest.n, forest.mst_a[keep], forest.mst_b[keep], forest.mst_w[keep], None
    )
    sub_leaves = np.flatnonzero(in_sub).astype(np.uint32)

    def fitch_on(f: Forest) -> int:
        tops = np.array([t for t in f.tops if t >= f.n or in_sub[t]], dtype=np.uint32)
        return int(schemes(f, leaves, tops)["fitch"].sum())

    return {
        "sequences": int(len(rows)),
        "storage_leaves": int(len(sub_leaves)),
        "nodes": n,
        "exact_candidates": len(ea),
        "exact_mst_weight": float(ed[exact["mst"]].sum()),
        "exact_components": len(exact["tops"]),
        "exact_s": exact_s,
        "seeded_candidates": len(ca),
        "seeded_largest_group": sstats["largest_group"],
        "seeded_mst_weight_exact": float(cd[seeded["mst"]].sum()),
        "seeded_components": len(seeded["tops"]),
        "seeded_s": seeded_s,
        "exact_mst_edges_in_candidates": len(exact_edges & cand_edges) / max(len(exact_edges), 1),
        "unreached": unreached,
        "storage_fitch_on_storage_tree": fitch_on(storage),
        "storage_fitch_on_sensitive_tree": fitch_on(sens),
    }


def measure(args: argparse.Namespace) -> None:
    """The sample's trees (``--members``, plus ``--extra`` representatives) -> ``--out``:

    - ``measure.json``: collapses, exact and staged spanning-tree weights, postings per
      scheme with and without floors and with floors on components only, cuts, chunking,
      density, one-topology-many-densities and the sensitive topology;
    - ``identity.tsv``, ``cuts.tsv``, ``changes.tsv``, ``gains.tsv``,
      ``prefix_postings.tsv``; ``kmers.parquet`` (per k-mer postings).
    """
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    members = load_members(args.members)
    if args.extra:
        extra = load_members(args.extra)
        members = (
            pl.concat([members, extra])
            .unique("protein_id", keep="first")
            .sort("cluster_rep", "protein_id")
        )
    seqs = seq_bytes(members)
    result: dict[str, Any] = {"params": {k: v for k, v in vars(args).items() if k != "func"}}
    cols = {}
    for tag, floors in (("floors", True), ("no_floors", False)):
        t0 = time.perf_counter()
        result[tag], cols[tag] = measure_variant(
            members, seqs, rates_of(args, floors), out, tag, args
        )
        result[tag]["wall_s"] = time.perf_counter() - t0
        print(tag, json.dumps(result[tag]["postings"]), flush=True)
    result["density"] = density(seqs, stored_sets(seqs, rates_of(args)), args.k)
    rows = subsample(members, args.dense_seqs)
    result["one_topology"] = one_topology(members, seqs, out, rows, args)
    col = cols["floors"]
    rows = subsample(members, args.sensitive_seqs)
    result["sensitive"] = sensitive_topology(seqs, out, col.leaf_of, rows, col.leaves, args)
    (out / "measure.json").write_text(json.dumps(result, indent=2, default=str) + "\n")


# ---------------------------------------------------------------------------- report


def report(args: argparse.Namespace) -> None:
    """The phase-0 gate: postings <= ``--max-postings`` x today's tier 1 and leaves <=
    ``--max-nodes`` x units, at matched t; plus the diff-scheme decision rule."""
    m = json.loads(Path(args.measure).read_text())
    tier1 = json.loads(Path(args.tier1).read_text())["stats"]
    p = m["floors"]["postings"]
    leaves, units = m["floors"]["collapse"]["leaves"], m["floors"]["units"]
    nodes = leaves + p["internal_nodes"]
    gate_postings = p["fitch"] <= args.max_postings * tier1["postings"]
    gate_nodes = leaves <= args.max_nodes * units
    rows = [
        ("gate", gate_postings and gate_nodes),
        ("gate_postings", gate_postings),
        ("gate_nodes", gate_nodes),
        ("tree_postings_fitch", p["fitch"]),
        ("tier1_postings", tier1["postings"]),
        ("postings_ratio", p["fitch"] / max(tier1["postings"], 1)),
        (
            "postings_ratio_no_floors",
            m["no_floors"]["postings"]["fitch"] / max(tier1["postings"], 1),
        ),
        (
            "postings_ratio_floors_on_components_only",
            m["floors"]["postings_floors_on_components_only"] / max(tier1["postings"], 1),
        ),
        ("leaves", leaves),
        ("internal_nodes", p["internal_nodes"]),
        ("units", units),
        ("leaves_per_unit", leaves / max(units, 1)),
        ("nodes_per_unit", nodes / max(units, 1)),
        ("sequence_mst_over_fitch", p["sequence_mst"] / max(p["fitch"], 1)),
        ("majority_over_fitch", p["majority"] / max(p["fitch"], 1)),
        ("add_only_over_fitch", p["add_only"] / max(p["fitch"], 1)),
        (
            "scheme_decision",
            "sequence_mst (if clearly cheaper to build: compare sequence_mst_s, fitch_s)"
            if p["sequence_mst"] <= 1.1 * p["fitch"]
            else "fitch_on_merge_tree",
        ),
    ]
    pl.DataFrame({"measure": [r[0] for r in rows], "value": [str(r[1]) for r in rows]}).write_csv(
        args.out, separator="\t"
    )


# ---------------------------------------------------------------------------- streaming


def batches(path: str, columns: list[str], rows: int = BATCH) -> Iterator[pl.DataFrame]:
    yield from pl.scan_parquet(path).select(columns).collect_batches(chunk_size=rows)


def sample(args: argparse.Namespace) -> None:
    """One bucket's sample clusters and their stored k-mers."""
    members = pl.read_parquet(args.members).filter(pl.col("cluster_rep") % args.sample == 0)
    members.write_parquet(f"{args.prefix}.members.parquet")
    sets = stored_sets(seq_bytes(load_members(f"{args.prefix}.members.parquet")), rates_of(args))
    np.save(f"{args.prefix}.hashes.npy", np.unique(sets.hashes))


def concat(args: argparse.Namespace) -> None:
    table = pl.concat([pl.read_parquet(p) for p in args.members])
    if args.max_rows and table.height > args.max_rows:  # the representatives sharing most
        table = table.sort("shared", "protein_id", descending=[True, False]).head(args.max_rows)
    table.sort("cluster_rep", "protein_id").write_parquet(args.out)
    if args.hashes:
        np.save(args.out_hashes, np.unique(np.concatenate([np.load(p) for p in args.hashes])))


def stream(args: argparse.Namespace) -> None:
    """One bucket's whole-release counts, pre-aggregated in the bucket and split by hash
    range (``--ranges``): ``digest`` (equal stored sets), ``signature`` (pigeonhole
    signatures, 1 in ``--sig-every`` kept by value), ``slice<rate>`` and ``dense<rate>``
    (k-mer holders in a hash window of width ``--slice`` at each rate, stored or all k-mers);
    plus the representatives outside the sample sharing stored k-mers with it."""
    rates = rates_of(args)
    core = np.load(args.core_hashes) if args.core_hashes else np.zeros(0, U64)
    parts: dict[str, list[pl.DataFrame]] = {}
    totals = {"sequences": 0, "residues": 0, "stored": 0, "empty": 0}
    extra = []
    width = args.slice
    for batch in batches(args.members, ["protein_id", "cluster_rep", "full_length", "sequence"]):
        batch = batch.with_columns(mask_adapters(pl.col("sequence")))
        seqs = seq_bytes(batch)
        sets = stored_sets(seqs, rates)
        totals["sequences"] += sets.n
        totals["residues"] += sum(len(s) for s in seqs)
        totals["stored"] += len(sets.hashes)
        totals["empty"] += int((sets.size == 0).sum())
        dig = digests(sets)[sets.size > 0, 0]
        parts.setdefault("digest", []).append(pl.DataFrame({"key": dig}))
        sig = signatures(sets, args.d + 1)
        key = _mix(sig["d0"].to_numpy() ^ sig["part"].to_numpy().astype(U64))
        key = key[key % U64(args.sig_every) == 0]
        parts.setdefault("signature", []).append(pl.DataFrame({"key": key}))
        for rate in args.slices:
            top = max(rate, width)  # the window (top - width, top]
            hi = _core.max_hash(top)
            lo = _core.max_hash(top - width) + 1 if top > width else 0
            seq, h = hashed(seqs, args.k, hi, min_hash=lo)
            parts.setdefault(f"dense{rate}", []).append(pl.DataFrame({"key": h}))
            stored = h <= sets.max_hash[seq]
            parts.setdefault(f"slice{rate}", []).append(pl.DataFrame({"key": h[stored]}))
        if len(core):
            is_rep = (batch["protein_id"] == batch["cluster_rep"]).to_numpy()
            outside = (batch["cluster_rep"] % args.sample != 0).to_numpy()
            hits = np.isin(sets.hashes, core)
            shared = np.bincount(sets.owner()[hits], minlength=sets.n)
            pick = np.flatnonzero(is_rep & outside & (shared >= args.min_shared))
            extra.append(batch[pick].with_columns(shared=pl.Series(shared[pick], dtype=pl.UInt32)))
    for kind, frames in parts.items():
        agg = pl.concat(frames).group_by("key").agg(count=pl.len().cast(pl.UInt32))
        agg = agg.with_columns(
            range=(pl.col("key") % args.ranges).cast(pl.UInt32), kind=pl.lit(kind)
        )
        parts[kind] = [agg]
    table = pl.concat([f[0] for f in parts.values()])
    for r, part in table.partition_by("range", as_dict=True).items():
        part.drop("range").write_parquet(f"{args.prefix}.r{r[0]}.parquet")
    for r in range(args.ranges):  # every range gets a file, possibly empty
        path = Path(f"{args.prefix}.r{r}.parquet")
        if not path.exists():
            table.clear().drop("range").write_parquet(path)
    if args.core_hashes:
        cols = {
            "protein_id": pl.Int64,
            "cluster_rep": pl.Int64,
            "full_length": pl.Boolean,
            "sequence": pl.String,
            "shared": pl.UInt32,
        }
        (pl.concat(extra) if extra else pl.DataFrame(schema=cols)).write_parquet(
            f"{args.prefix}.extra.parquet"
        )
    Path(f"{args.prefix}.totals.json").write_text(json.dumps(totals) + "\n")


def reduce(args: argparse.Namespace) -> None:
    """One hash range of every bucket -> per kind, histogram of group sizes and the largest."""
    table = pl.read_parquet(args.parts).group_by("kind", "key").agg(pl.sum("count"))
    hist = table.group_by("kind", size="count").len("n")
    top = table.sort("count", descending=True).group_by("kind").head(args.top)
    hist.write_parquet(f"{args.prefix}.hist.parquet")
    top.write_parquet(f"{args.prefix}.top.parquet")


def combine(args: argparse.Namespace) -> None:
    hist = pl.read_parquet(args.hist).group_by("kind", "size").agg(pl.sum("n")).sort("kind", "size")
    hist.write_csv("stream.tsv", separator="\t")
    totals: dict[str, int] = {}
    for path in args.totals:
        for k, v in json.loads(Path(path).read_text()).items():
            totals[k] = totals.get(k, 0) + v
    summary: dict[str, Any] = {"totals": totals}
    for (kind,), part in hist.partition_by("kind", as_dict=True).items():
        groups, members = int(part["n"].sum()), int((part["size"] * part["n"]).sum())
        summary[str(kind)] = {
            "groups": groups,
            "members": members,
            "largest": int(part["size"].max()),  # type: ignore[arg-type]
            "distinct_over_total": groups / max(members, 1),
        }
    nonempty = totals["sequences"] - totals["empty"]
    summary["equal_set_collapse_share"] = 1 - summary["digest"]["groups"] / max(nonempty, 1)
    summary["largest_groups"] = (
        pl.read_parquet(args.tops)
        .sort("count", descending=True)
        .group_by("kind")
        .head(args.top)
        .to_dicts()
    )
    Path("stream.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")


# ---------------------------------------------------------------------------- fixture


def fixture(args: argparse.Namespace) -> None:
    """Synthetic families: a root protein per family; members at 100-80% identity, a 90%
    cluster around the root (the others singletons), and fragments of members as partial
    sequences; written as ``--buckets`` members tables split by cluster."""
    rng = np.random.default_rng(args.seed)
    aa = np.frombuffer(b"ACDEFGHIKLMNPQRSTVWY", dtype=np.uint8)
    rows = []
    pid = 1
    for _ in range(args.families):
        root = rng.choice(aa, rng.integers(120, 600))
        rep = pid
        for ident in rng.choice([1.0, 0.995, 0.99, 0.97, 0.95, 0.92, 0.85, 0.8], args.members):
            seq = root.copy()
            sites = rng.random(len(seq)) > ident
            seq[sites] = rng.choice(aa, sites.sum())
            cluster = rep if ident >= 0.9 else pid
            rows.append((pid, cluster, True, seq.tobytes().decode()))
            pid += 1
            if rng.random() < 0.5:  # a fragment of it
                start = rng.integers(0, len(seq) // 2)
                frag = seq[start : start + rng.integers(40, len(seq) - start + 1)]
                rows.append((pid, cluster, False, frag.tobytes().decode()))
                pid += 1
    table = pl.DataFrame(
        rows, schema=["protein_id", "cluster_rep", "full_length", "sequence"], orient="row"
    ).with_columns(bucket=pl.col("cluster_rep").hash(1) % args.buckets)
    for b in range(args.buckets):
        table.filter(pl.col("bucket") == b).drop("bucket").sort(
            "cluster_rep", "protein_id"
        ).write_parquet(f"members.bucket{b}.parquet")


# ---------------------------------------------------------------------------- cli


def floats(text: str) -> list[float]:
    return [float(x) for x in text.split(",")]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(required=True)

    def command(name: str, func: Any) -> argparse.ArgumentParser:
        p = sub.add_parser(name)
        p.set_defaults(func=func)
        p.add_argument("--k", type=int, default=11)
        p.add_argument("--t-base", type=float, default=0.02)
        p.add_argument("--n-min", type=int, default=16)
        p.add_argument("--t-cap", type=float, default=0.2)
        p.add_argument("--levels", type=int, default=8)
        p.add_argument("--d", type=int, default=3, help="pigeonhole stage: weights <= d")
        p.add_argument("--seed", type=int, default=1)
        return p

    p = command("fixture", fixture)
    p.add_argument("--families", type=int, default=60)
    p.add_argument("--members", type=int, default=12)
    p.add_argument("--buckets", type=int, default=2)

    p = command("sample", sample)
    p.add_argument("--members", required=True)
    p.add_argument("--sample", type=int, required=True)
    p.add_argument("--prefix", required=True)

    p = command("concat", concat)
    p.add_argument("--members", nargs="+", required=True)
    p.add_argument("--hashes", nargs="*", default=[])
    p.add_argument("--out", required=True)
    p.add_argument("--out-hashes", default="core_hashes.npy")
    p.add_argument("--max-rows", type=int, default=0, help="keep the rows sharing most")

    p = command("stream", stream)
    p.add_argument("--members", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--core-hashes", default="")
    p.add_argument("--sample", type=int, default=1)
    p.add_argument("--min-shared", type=int, default=4)
    p.add_argument("--ranges", type=int, default=16)
    p.add_argument("--sig-every", type=int, default=1)
    p.add_argument("--slice", type=float, default=1 / 1024)
    p.add_argument("--slices", type=floats, default=[0.0, 0.06, 0.2])

    p = command("reduce", reduce)
    p.add_argument("--parts", nargs="+", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--top", type=int, default=20)

    p = command("combine", combine)
    p.add_argument("--hist", nargs="+", required=True)
    p.add_argument("--tops", nargs="+", required=True)
    p.add_argument("--top", type=int, default=20)
    p.add_argument("--totals", nargs="+", required=True)

    p = command("measure", measure)
    p.add_argument("--members", nargs="+", required=True)
    p.add_argument("--extra", nargs="*", default=[])
    p.add_argument("--out", required=True)
    p.add_argument(
        "--group-cap", type=int, default=64, help="stage 3: all pairs up to this group size"
    )
    p.add_argument(
        "--posting-cap", type=int, default=1000, help="stage 4: longer prefix postings fall back"
    )
    p.add_argument("--sn-orders", type=int, default=4)
    p.add_argument("--sn-window", type=int, default=2)
    p.add_argument("--max-groups", type=int, default=16, help="today's promiscuity cut")
    p.add_argument("--identity-pairs", type=int, default=20000)
    p.add_argument("--chunk-sizes", type=floats, default=[1e3, 1e4, 1e5, 1e6])
    p.add_argument("--dense-seqs", type=int, default=50000)
    p.add_argument("--sensitive-seqs", type=int, default=50000)
    p.add_argument("--seed-alphabet", default="murphy10")
    p.add_argument("--seed-k", type=int, default=10)
    p.add_argument("--seeds", type=int, default=20)
    p.add_argument("--sketch-rate", type=float, default=0.25)

    p = command("report", report)
    p.add_argument("--measure", required=True)
    p.add_argument("--tier1", required=True, help="today's index meta.json on the same proteins")
    p.add_argument("--max-postings", type=float, default=1.5)
    p.add_argument("--max-nodes", type=float, default=2.0)
    p.add_argument("--out", default="gate.tsv")

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
