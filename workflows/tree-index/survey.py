"""Tree index, floor survey (docs/tree-index-plan.md, Progress log: forest floors, pieces by
coherence, floors on every piece, size options).

One topology per sample: phase 0's exact storage tree on per-sequence floored sets (the full
hierarchy, kept whole). On it, for each configuration (t\\_base, n\\_min, tau, n\\_single, n\\_1):

- **base band:** every leaf's k-mers with hash <= t\\_base, diffed by Fitch over the whole
  tree;
- **piece floors:** top-down from each component, a node whose floor is at least ``tau``
  coherent is a piece (``tau`` 0: the components). Its floor is tier 1's, t\\_p = max(t\\_base,
  min(t\\_cap, n\\_min / its union's k-mers)), the n\\_min candidates above t\\_base held by most
  members, posted once at its root. Coherence is the mean share of the floor a member holds.
  Every piece is floored, single leaves included;
- **single-leaf pieces:** a floor of their own ``n_single`` lowest hashes above t\_base (no one
  to share it with, so each is a posting); ``n_single`` = n\_min is the uniform rule;
- **leaf minimum:** each leaf topped up to ``n_1`` stored k-mers of its own (base, floor
  k-mers it holds, then its lowest other hashes), each top-up k-mer one posting.

Writes ``survey.tsv`` (one row per configuration, with tier-1 ratios and release-scale
projections) and ``survey.json`` (hierarchy reach of proteins without own base k-mers, tier
1's own floor coherence). A separate script, so phase0.py and the cached tasks stay as they are.
"""

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np
import phase0 as p0
import polars as pl
from numpy.typing import NDArray

from kmer_functional_profiler import _core

U64 = np.uint64


class Tree:
    """The storage tree and every protein's k-mers above the lowest t\\_base, in leaf order."""

    def __init__(self, members: pl.DataFrame, args: argparse.Namespace) -> None:
        self.args = args
        seqs = p0.seq_bytes(members)
        self.n_proteins = len(seqs)
        sets = p0.stored_sets(
            seqs,
            p0.Rates(args.k, args.topology_t_base, args.topology_n_min, args.t_cap, args.levels),
        )
        col, pairs = p0.collapse(sets, members["full_length"].to_numpy())
        self.col, self.leaves = col, col.leaves
        t0 = time.perf_counter()
        self.forest = p0.spanning_forest(
            self.leaves.n, pairs["a"], pairs["b"], pairs["weight"], self.leaves.size
        )
        self.build_s = time.perf_counter() - t0
        f = self.forest
        self.lay = _core.tree_layout(f.n, f.left, f.right, f.tops)
        pos_of_leaf = np.empty(f.n, dtype=np.int64)
        pos_of_leaf[self.lay["order"]] = np.arange(f.n)
        # proteins with a leaf, in leaf order: one representative per leaf (its longest)
        lengths = np.array([len(s) for s in seqs])
        has = np.flatnonzero(col.leaf_of >= 0)
        rep = (
            pl.DataFrame({"seq": has, "leaf": col.leaf_of[has], "len": lengths[has]})
            .sort("leaf", "len", "seq", descending=[False, True, False])
            .group_by("leaf", maintain_order=True)
            .first()
            .sort("leaf")["seq"]
            .to_numpy()
        )
        self.rep_by_pos = rep[self.lay["order"]]  # the leaf at each position's representative
        self.members_per_leaf = np.bincount(col.leaf_of[has], minlength=f.n)
        self.leaf_of = col.leaf_of
        self.n_kmers = _core.distinct_kmers(seqs, np.arange(len(seqs), dtype=np.uint32), args.k)
        self.n_valid = np.zeros(len(seqs), dtype=np.int64)
        self.n_valid[self.n_kmers["group"]] = self.n_kmers["n_kmers"]
        # every representative's k-mers up to t_cap, rows grouped by leaf position
        reps = [seqs[i] for i in self.rep_by_pos]
        seq, h = p0.hashed(reps, args.k, _core.max_hash(args.t_cap))
        self.pos, self.h = seq, h  # sorted by (position, hash)
        self.row = np.searchsorted(seq, np.arange(f.n + 1))
        self.seqs = seqs


def base_postings(tree: Tree, t_base: float) -> tuple[int, int, NDArray[np.int64]]:
    """Fitch postings and distinct k-mers of the base band over the whole tree; and each
    leaf position's base size."""
    mh = U64(_core.max_hash(t_base))
    keep = tree.h <= mh
    pos, h = tree.pos[keep], tree.h[keep]
    n = tree.forest.n
    leaf = tree.lay["order"][pos]  # leaf id of each row
    o = np.lexsort((h, leaf))
    leaf, h = leaf[o], h[o]
    sets = p0.Sets(p0.csr(leaf, h, n), h, np.zeros(n, np.uint8), np.array([mh], dtype=U64))
    fitch = int(p0.schemes(tree.forest, sets)["fitch"].sum())
    size_by_pos = np.bincount(pos, minlength=n)
    return fitch, len(np.unique(h)), size_by_pos


def pieces(
    tree: Tree, t_base: float, n_min: int, tau: float
) -> tuple[list[tuple[int, int, NDArray[np.uint64], float]], int]:
    """Coherence pieces: (lo, hi, floor hashes, coherence) per piece, and nodes evaluated."""
    a = tree.args
    lo_base = U64(_core.max_hash(t_base))
    f, lay = tree.forest, tree.lay
    out, evaluated = [], 0
    stack = [int(t) for t in f.tops]
    while stack:
        x = stack.pop()
        lo, hi = int(lay["lo"][x]), int(lay["hi"][x])
        hs = tree.h[tree.row[lo] : tree.row[hi]]
        hs = hs[hs > lo_base]
        uniq, counts = np.unique(hs, return_counts=True)
        evaluated += 1
        n_kmers = max(len(uniq) / (a.t_cap - t_base), 1.0)
        t = max(t_base, min(a.t_cap, n_min / n_kmers))
        sel = uniq <= U64(_core.max_hash(t))
        best = np.argsort(-counts[sel], kind="stable")[:n_min]
        floor, held = uniq[sel][best], counts[sel][best]
        leaves = hi - lo
        coh = float(held.mean() / leaves) if len(held) else 1.0
        if x < f.n or tau <= 0 or coh >= tau:
            out.append((lo, hi, floor, coh))
            continue
        i = x - f.n
        stack += [int(f.left[i]), int(f.right[i])]
    return out, evaluated


def configuration(
    tree: Tree,
    base: tuple[int, int, NDArray[np.int64]],
    t_base: float,
    n_min: int,
    tau: float,
    n_singles: list[int],
    n_1s: list[int],
) -> list[dict[str, Any]]:
    """Rows for one piece set: each single-leaf floor size (``n_single``, at most n_min) and
    leaf minimum (``n_1``)."""
    fitch, base_keys, base_size = base
    t0 = time.perf_counter()
    got, evaluated = pieces(tree, t_base, n_min, tau)
    lo_base = U64(_core.max_hash(t_base))
    n = tree.forest.n
    held = np.zeros(n, dtype=np.int64)  # floor k-mers each leaf position holds
    avail = np.zeros(n, dtype=np.int64)  # its other hashes above t_base (up to t_cap)
    multi_postings, multi_keys, cohs = 0, [], []
    single_pos, single_len, single_keys = [], [], []
    for lo, hi, floor, coh in got:
        r0, r1 = tree.row[lo], tree.row[hi]
        hs, pos = tree.h[r0:r1], tree.pos[r0:r1]
        above = hs > lo_base
        if hi - lo == 1:  # its floor: its own lowest hashes above t_base, best first
            single_pos.append(lo)
            single_len.append(len(floor))
            single_keys.append(floor)
            avail[lo] = int(above.sum())
            continue
        multi_postings += len(floor)
        multi_keys.append(floor)
        cohs.append(coh)
        hit = above & np.isin(hs, floor)
        held[lo:hi] = np.bincount(pos[hit] - lo, minlength=hi - lo)
        avail[lo:hi] = np.bincount(pos[above & ~hit] - lo, minlength=hi - lo)
    pieces_s = time.perf_counter() - t0
    sp, sl = np.array(single_pos, dtype=np.int64), np.array(single_len, dtype=np.int64)
    rep_valid = tree.n_valid[tree.rep_by_pos]
    rows = []
    for n_single in sorted({min(x, n_min) for x in n_singles}):
        kept = np.minimum(sl, n_single)
        held_s, avail_s = held.copy(), avail.copy()
        held_s[sp] = kept
        avail_s[sp] -= kept
        floor_postings = multi_postings + int(kept.sum())
        keys = [*multi_keys, *(k[:n_single] for k in single_keys)]
        n_keys = len(np.unique(np.concatenate(keys))) if keys else 0
        own = base_size + held_s
        for n_1 in n_1s:
            need = np.maximum(0, n_1 - own)
            # hashes above t_cap are available too: a protein's valid k-mers cap the top-up
            topup = np.minimum(need, np.maximum(rep_valid - own, 0))
            covered = (own + topup) > 0
            # proteins with no k-mer up to t_cap have no leaf: under a minimum each becomes a
            # leaf of its own, of its n_1 lowest hashes (nothing to share them with)
            leafless = tree.n_valid[tree.leaf_of < 0]
            new_leaves = np.minimum(leafless, n_1)
            topup_total = int(topup.sum()) + int(new_leaves.sum())
            blind = int(tree.members_per_leaf[tree.lay["order"][~covered]].sum()) + int(
                (new_leaves == 0).sum()
            )
            rows.append(
                {
                    "t_base": t_base,
                    "n_min": n_min,
                    "tau": tau,
                    "n_single": n_single,
                    "n_1": n_1,
                    "pieces": len(got),
                    "single_leaf_pieces": len(sp),
                    "pieces_mean_coherence": float(np.mean(cohs)) if cohs else None,
                    "base_postings": fitch,
                    "floor_postings": floor_postings,
                    "topup_postings": topup_total,
                    "leafless_proteins": len(leafless),
                    "topup_beyond_t_cap": int(np.maximum(topup - avail_s, 0).sum()),
                    "postings": fitch + floor_postings + topup_total,
                    "keys_base": base_keys,
                    "keys_floor": n_keys,
                    "blind_proteins": blind,
                    "nodes_evaluated": evaluated,
                    "pieces_s": pieces_s,
                }
            )
    return rows


def reach(tree: Tree, t_base: float) -> dict[str, int]:
    """How proteins without an own base k-mer can reach tier 2 through the hierarchy."""
    f, lay = tree.forest, tree.lay
    comp = np.zeros(f.n, dtype=np.int64)
    for t in f.tops:
        lo, hi = lay["lo"][t], lay["hi"][t]
        comp[lay["order"][lo:hi]] = hi - lo
    dense = p0.dense_sets(tree.seqs, tree.args.k, t_base).size
    has = tree.leaf_of >= 0
    blind = dense == 0
    alone = np.zeros(len(blind), dtype=bool)
    alone[has] = comp[tree.leaf_of[has]] == 1
    shared = np.zeros(len(blind), dtype=bool)
    shared[has] = tree.members_per_leaf[tree.leaf_of[has]] > 1
    return {
        "no_own_base_kmer": int(blind.sum()),
        "in_component_with_relatives": int((blind & has & ~alone).sum()),
        "alone_sharing_a_leaf": int((blind & alone & shared).sum()),
        "alone": int((blind & alone & ~shared).sum()),
        "no_leaf": int((blind & ~has).sum()),
        "no_valid_kmer": int((tree.n_valid == 0).sum()),
    }


def tier1_coherence(
    members: pl.DataFrame, tree: Tree, t_base: float, n_min: int
) -> dict[str, float]:
    """Tier 1's floor coherence by the same measure, its units being the clusters."""
    seqs = tree.seqs
    lo_base = U64(_core.max_hash(tree.args.t_cap))
    seq, h = p0.hashed(seqs, tree.args.k, int(lo_base), min_hash=_core.max_hash(t_base) + 1)
    units = members["cluster_rep"].to_numpy()[seq]
    table = pl.DataFrame({"unit": units, "seq": seq, "hash": h})
    size = pl.DataFrame({"unit": members["cluster_rep"]}).group_by("unit").len("members")
    n = table.group_by("unit").agg(uniq=pl.col("hash").n_unique())
    t = n.with_columns(
        t=pl.max_horizontal(
            pl.lit(t_base),
            pl.min_horizontal(
                pl.lit(tree.args.t_cap),
                n_min / (pl.col("uniq") / (tree.args.t_cap - t_base)).clip(1),
            ),
        )
    ).with_columns(mh=pl.col("t").map_elements(_core.max_hash, return_dtype=pl.UInt64))
    coh = (
        table.join(t.select("unit", "mh"), on="unit")
        .filter(pl.col("hash") <= pl.col("mh"))
        .group_by("unit", "hash")
        .agg(holders=pl.col("seq").n_unique())
        .sort(["unit", "holders", "hash"], descending=[False, True, False])
        .filter(pl.int_range(pl.len()).over("unit") < n_min)
        .group_by("unit")
        .agg(pl.col("holders").mean())
        .join(size, on="unit")
        .filter(pl.col("members") > 1)
        .with_columns(c=pl.col("holders") / pl.col("members"))
    )
    return {
        "units": coh.height,
        "mean": float(coh["c"].mean() or 0),
        "median": float(coh["c"].median() or 0),
    }  # type: ignore[arg-type]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--members", nargs="+", required=True)
    parser.add_argument("--extra", nargs="*", default=[])
    parser.add_argument("--out", required=True)
    parser.add_argument("--tier1", nargs="*", default=[], help="tier-1 meta.json files to compare")
    parser.add_argument("--stream", default="", help="stream.json, for release-scale projections")
    parser.add_argument("--k", type=int, default=11)
    parser.add_argument("--t-cap", type=float, default=0.2)
    parser.add_argument("--levels", type=int, default=8)
    parser.add_argument("--topology-t-base", type=float, default=0.02)
    parser.add_argument("--topology-n-min", type=int, default=16)
    parser.add_argument("--t-bases", type=p0.floats, default=[0.01, 0.02])
    parser.add_argument("--n-mins", type=lambda s: [int(x) for x in s.split(",")], default=[8, 16])
    parser.add_argument("--taus", type=p0.floats, default=[0.0, 0.3, 0.6, 0.9])
    parser.add_argument(
        "--n-singles", type=lambda s: [int(x) for x in s.split(",")], default=[0, 2, 4, 16]
    )
    parser.add_argument(
        "--n-1s", type=lambda s: [int(x) for x in s.split(",")], default=[0, 1, 2, 4]
    )
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    members = p0.load_members(args.members)
    if args.extra:
        extra = p0.load_members(args.extra)
        members = (
            pl.concat([members, extra])
            .unique("protein_id", keep="first")
            .sort("cluster_rep", "protein_id")
        )
    tree = Tree(members, args)
    rows: list[dict[str, Any]] = []
    summary: dict[str, Any] = {
        "proteins": tree.n_proteins,
        "units": int(members["cluster_rep"].n_unique()),
        "leaves": tree.forest.n,
        "internal_nodes": len(tree.forest.left),
        "components": len(tree.forest.tops),
        "topology_build_s": tree.build_s,
        "reach": {},
        "tier1_coherence": {},
    }
    for t_base in args.t_bases:
        base = base_postings(tree, t_base)
        summary["reach"][str(t_base)] = reach(tree, t_base)
        for n_min in args.n_mins:
            summary["tier1_coherence"][f"{t_base}/{n_min}"] = tier1_coherence(
                members, tree, t_base, n_min
            )
            for tau in args.taus:
                for r in configuration(tree, base, t_base, n_min, tau, args.n_singles, args.n_1s):
                    rows.append(r)
                    print(json.dumps(r), flush=True)
    table = pl.DataFrame(rows)
    # tier 1 at matching rates (and the standard tier for every row)
    tiers = []
    for path in args.tier1:
        meta = json.loads(Path(path).read_text())
        prm = meta["params"]
        tiers.append(
            {"t_base": prm["t_base"], "n_min": prm["n_min"], "postings": meta["stats"]["postings"]}
        )
    standard = next((t for t in tiers if t["t_base"] == 0.02 and t["n_min"] == 16), None)
    if standard:
        table = table.with_columns(over_standard=pl.col("postings") / standard["postings"])
    if tiers:
        matched = pl.DataFrame(tiers).rename({"postings": "tier1_matched"})
        table = table.join(matched, on=["t_base", "n_min"], how="left").with_columns(
            over_matched=pl.col("postings") / pl.col("tier1_matched")
        )
    if args.stream:
        release = json.loads(Path(args.stream).read_text())["totals"]["sequences"]
        scale = release / tree.n_proteins
        summary["release_scale"] = scale
        table = table.with_columns(
            release_postings=pl.col("postings") * scale,
            release_keys=(pl.col("keys_base") + pl.col("keys_floor")) * scale,
            release_pieces=pl.col("pieces") * scale,
        )
        summary["release_leaves"] = tree.forest.n * scale
    summary["tier1"] = tiers
    table.write_csv(out / "survey.tsv", separator="\t")
    (out / "survey.json").write_text(json.dumps(summary, indent=1) + "\n")


if __name__ == "__main__":
    main()
