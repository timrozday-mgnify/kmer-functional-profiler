"""Tree index: split proteins into parts for a parallel build (``build.nf``).

Proteins are linked when they can share a storage-tree edge at the topology rates (k 11 stored
sets, ``_core.overlap_components``: the pairs are unioned as found, not kept), or have the same
sequence (so equal dense sets, which share a
leaf, never straddle parts). Connected components are packed, largest first, into parts of
about ``--part-proteins`` proteins; a component larger than that is a part of its own. Each part
is built alone (``build.py``) and the parts are merged (``merge.py``). Edges the dense topology
would find between components are lost, so a k-mer shared across parts is posted in each.

Writes ``part0000.parquet``... (the members table, per part) and ``partition.json``.
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


def components(members: pl.DataFrame, args: argparse.Namespace) -> tuple[np.ndarray, dict]:
    """Each protein's component, and stats."""
    t0 = time.perf_counter()
    seqs = p0.seq_bytes(members)
    rates = p0.Rates(args.k, args.topology_t_base, args.topology_n_min, args.t_cap, args.levels)
    sets = p0.stored_sets(seqs, rates)
    linked = _core.overlap_components(*sets.kernel()).astype(np.int64)
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
    sizes = np.bincount(comp)
    return comp, {
        "proteins": n,
        "same_sequence": twins.height,
        "components": int(n_comp),
        "largest_component": int(sizes.max()) if n else 0,
        "components_s": time.perf_counter() - t0,
    }


def pack(comp: np.ndarray, target: int) -> np.ndarray:
    """Each protein's part: components largest first, cut every ``target`` proteins."""
    sizes = np.bincount(comp)
    order = np.argsort(-sizes, kind="stable")
    start = np.cumsum(sizes[order]) - sizes[order]  # proteins before each component
    part_of = np.empty(len(sizes), dtype=np.int64)
    part_of[order] = start // target
    # renumber densely (a large component can skip numbers)
    _, part_of = np.unique(part_of, return_inverse=True)
    return part_of[comp]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--members", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--part-proteins", type=int, default=2_000_000)
    parser.add_argument("--k", type=int, default=11)
    parser.add_argument("--t-cap", type=float, default=0.2)
    parser.add_argument("--levels", type=int, default=8)
    parser.add_argument("--topology-t-base", type=float, default=0.02)
    parser.add_argument("--topology-n-min", type=int, default=16)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    members = p0.load_members(args.members)
    comp, stats = components(members, args)
    part = pack(comp, args.part_proteins)
    # the raw columns (adapters unmasked): build.py masks them again on load
    raw = pl.read_parquet(args.members, columns=["protein_id", "cluster_rep", "full_length",
                                                 "sequence"])  # fmt: skip
    raw = raw.sort("cluster_rep", "protein_id").with_columns(part=part)
    sizes = []
    for (p,), rows in raw.partition_by("part", as_dict=True, maintain_order=True).items():
        rows.drop("part").write_parquet(out / f"part{p:04d}.parquet")
        sizes.append(rows.height)
    stats |= {"parts": len(sizes), "part_proteins": sorted(sizes, reverse=True)}
    (out / "partition.json").write_text(json.dumps(stats, indent=1) + "\n")
    print(json.dumps(stats | {"part_proteins": sizes[:5]}))


if __name__ == "__main__":
    main()
