"""Tree index: split proteins into parts for a parallel build (``build.nf``).

Proteins are linked when they can share a storage-tree edge at the topology rates (k 11 stored
sets, ``_core.overlap_components``: the pairs are unioned as found, not kept), or have the same
sequence (so equal dense sets, which share a
leaf, never straddle parts). Connected components are packed, largest first, into parts of
about ``--part-proteins`` proteins. A component larger than that (single linkage chains families
together once clusters are dense: the 1-in-100 subset's largest outgrew BUILD's memory) is cut
between atoms, an MGnify90 cluster with any protein of an equal sequence, taken in order of each
atom's least stored hash, so atoms sharing k-mers tend to share a part. Each part is built alone
(``build.py``) and the parts are merged (``merge.py``). Edges the dense topology would find between
parts are lost, so a k-mer shared across parts is posted in each.

Writes ``part0000.parquet``... (the members table, per part) and ``partition.json``, each name
after ``--prefix``. At release, one call per members bucket (``build.nf --partition_per_file``):
buckets hold whole MGnify90 clusters, so only edges between clusters in different buckets are
lost, and their k-mers are posted once per bucket.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import phase0 as p0
import polars as pl
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from kmer_functional_profiler import _core


def components(
    members: pl.DataFrame, args: argparse.Namespace
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """Each protein's component, atom (its cluster and equal sequences) and least stored hash
    (``2^64 - 1`` for none), and stats."""
    t0 = time.perf_counter()
    seqs = p0.seq_bytes(members)
    rates = p0.Rates(args.k, args.topology_t_base, args.topology_n_min, args.t_cap, args.levels)
    sets = p0.stored_sets(seqs, rates)
    linked = _core.overlap_components(*sets.kernel()).astype(np.int64)
    key = np.full(sets.n, np.iinfo(np.uint64).max, dtype=np.uint64)
    held = sets.size > 0
    key[held] = sets.hashes[sets.offsets[:-1][held]]  # sets are sorted: the first is the least
    del seqs, sets
    same = members.select(h=pl.col("sequence").hash()).with_row_index("i")
    first = same.group_by("h").agg(pl.col("i").min().alias("j"))
    twins = same.join(first, on="h").filter(pl.col("i") != pl.col("j"))
    n = members.height
    # each protein to its component's first set, and equal sequences to each other
    a = np.concatenate([np.arange(n), twins["i"].to_numpy().astype(np.int64)])
    b = np.concatenate([linked, twins["j"].to_numpy().astype(np.int64)])
    graph = coo_matrix((np.ones(len(a), dtype=np.int8), (a, b)), shape=(n, n))
    n_comp, comp = connected_components(graph, directed=False)
    # atoms: each protein to its cluster (node n + c), and equal sequences to each other
    _, cluster = np.unique(members["cluster_rep"].to_numpy(), return_inverse=True)
    m = n + int(cluster.max()) + 1 if n else 0
    a = np.concatenate([np.arange(n), twins["i"].to_numpy().astype(np.int64)])
    b = np.concatenate([n + cluster, twins["j"].to_numpy().astype(np.int64)])
    _, atom = connected_components(
        coo_matrix((np.ones(len(a), dtype=np.int8), (a, b)), shape=(m, m)), directed=False
    )
    sizes = np.bincount(comp)
    return (
        comp,
        atom[:n],
        key,
        {
            "proteins": n,
            "same_sequence": twins.height,
            "components": int(n_comp),
            "largest_component": int(sizes.max()) if n else 0,
            "components_s": time.perf_counter() - t0,
        },
    )


def pack(comp: np.ndarray, target: int, atom: np.ndarray, key: np.ndarray) -> np.ndarray:
    """Each protein's part: components largest first, cut every ``target`` proteins; a
    component larger than ``target`` in pieces, its atoms in order of their least ``key``.
    ponytail: one least hash per atom (a MinHash of one) places related atoms together only
    by chance agreement; ordering by the component's own edges if parts duplicate too much."""
    sizes = np.bincount(comp)
    big = sizes[comp] > target
    # pieces: a whole component, or an atom of a large one
    _, piece = np.unique(np.where(big, len(sizes) + atom, comp), return_inverse=True)
    n_piece = int(piece.max()) + 1 if len(piece) else 0
    p_comp = np.empty(n_piece, dtype=np.int64)
    p_comp[piece] = comp
    p_key = np.full(n_piece, np.iinfo(np.uint64).max, dtype=np.uint64)
    np.minimum.at(p_key, piece, key)
    p_size = np.bincount(piece, minlength=n_piece)
    order = np.lexsort((p_key, p_comp, -sizes[p_comp]))
    start = np.cumsum(p_size[order]) - p_size[order]  # proteins before each piece
    part_of = np.empty(n_piece, dtype=np.int64)
    part_of[order] = start // target
    # renumber densely (a large piece can skip numbers)
    _, part_of = np.unique(part_of, return_inverse=True)
    return part_of[piece]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--members", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--part-proteins", type=int, default=2_000_000)
    parser.add_argument("--prefix", default="", help="before each output's name (a bucket's)")
    parser.add_argument("--k", type=int, default=11)
    parser.add_argument("--t-cap", type=float, default=0.2)
    parser.add_argument("--levels", type=int, default=8)
    parser.add_argument("--topology-t-base", type=float, default=0.02)
    parser.add_argument("--topology-n-min", type=int, default=16)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    members = p0.load_members(args.members)
    comp, atom, key, stats = components(members, args)
    part = pack(comp, args.part_proteins, atom, key)
    stats["split_components"] = int((np.bincount(comp) > args.part_proteins).sum())
    # the raw columns (adapters unmasked): build.py masks them again on load
    raw = pl.read_parquet(args.members, columns=["protein_id", "cluster_rep", "full_length",
                                                 "sequence"])  # fmt: skip
    raw = raw.sort("cluster_rep", "protein_id").with_columns(part=part)
    sizes = []
    for (p,), rows in raw.partition_by("part", as_dict=True, maintain_order=True).items():
        rows.drop("part").write_parquet(out / f"{args.prefix}part{p:04d}.parquet")
        sizes.append(rows.height)
    stats |= {"parts": len(sizes), "part_proteins": sorted(sizes, reverse=True)}
    (out / f"{args.prefix}partition.json").write_text(json.dumps(stats, indent=1) + "\n")
    print(json.dumps(stats | {"part_proteins": sizes[:5]}))


if __name__ == "__main__":
    main()
