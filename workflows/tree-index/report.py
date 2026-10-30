"""Tree index, phase 3: reporting, recall first (docs/tree-index-plan.md, Querying a metagenome,
step 3, and Recall-first candidates and neighbours), on ``query.py``'s counts.

- **AAI** per hit leaf and node: containment^(1/k), with a Clopper-Pearson ``level`` interval on
  its hit k-mers out of |S\\_v| (``aai``, ``aai_lo``, ``aai_hi``). Naive: depth below full
  coverage lowers containment as divergence does; phase 4 calibrates it against spiked truth.
- **Candidates:** every leaf with at least ``min_hits`` hit k-mers is a seed. Each adds its
  single-linkage cluster on the spanning tree at R\\_v = d(aai\\_lo) + d(1 - delta), where
  d(p) = 2 (len\\_v - k\\_t + 1)(1 - p^k\\_t) is the expected |A Δ B| of two proteins' dense
  k\\_t-mer sets (the topology's units, ``--leaves proteins``) at identity p. A cluster is walked
  nearest first in minimax distance (Prim from the seed), at most ``max_neighbours`` leaves
  beyond the seed. Candidate positions are written as ranges, and their proteins as members,
  for round 2 and the pileup's reference.
- **Gather** over the seeds: take the leaf with the most remaining hit k-mers, scaled by its
  windows over |S\\_v| (estimated shared windows, as today's gather scales by 1 / t\\_g). Leaves
  in its component sharing a remaining k-mer with it, whose score is within ``tie_z`` standard
  deviations of its (Poisson on its count), are tied, and the group is reported as its smallest
  common subtree. The group's k-mers are removed; repeat while the best has ``min_hits``.

Writes to ``--out``: ``leaves.parquet`` and ``nodes.parquet`` (the query's, with AAI);
``candidates.parquet`` (lo, hi: position ranges); ``members.parquet`` (protein\\_id,
cluster\\_rep, pos, seed); ``gather.parquet`` (rank, pos, score, kmers, tied, node, lo, hi, aai,
aai\\_lo, aai\\_hi; node -1 when the group is one leaf); and ``report.json``.
"""

import argparse
import heapq
import json
import time
from pathlib import Path

import numpy as np
import polars as pl
import query
from numpy.typing import NDArray
from scipy.stats import beta


def with_aai(t: pl.DataFrame, k: int, level: float) -> pl.DataFrame:
    """containment^(1/k), with a Clopper-Pearson interval on hit\\_kmers of size."""
    h = t["hit_kmers"].to_numpy().astype(np.float64)
    n = np.maximum(t["size"].to_numpy().astype(np.float64), h)
    q = (1 - level) / 2
    with np.errstate(all="ignore"):
        lo = np.where(h > 0, beta.ppf(q, h, n - h + 1), 0.0)
        hi = np.where(h < n, beta.ppf(1 - q, h + 1, n - h), 1.0)
    root = lambda c: np.nan_to_num(c, nan=0.0) ** (1 / k)  # noqa: E731
    return t.with_columns(
        aai=pl.Series(root(h / np.maximum(n, 1))),
        aai_lo=pl.Series(root(lo)),
        aai_hi=pl.Series(root(hi)),
    )


class Spanning:
    """The spanning tree over positions, as adjacency lists."""

    def __init__(self, path: Path, n_pos: int) -> None:
        a, b, w = (np.load(path / f"edge_{x}.npy").astype(np.int64) for x in "abw")
        src = np.concatenate([a, b])
        o = np.argsort(src, kind="stable")
        self.dst = np.concatenate([b, a])[o]
        self.w = np.concatenate([w, w])[o]
        self.start = np.searchsorted(src[o], np.arange(n_pos + 1))

    def cluster(self, v: int, radius: float, cap: int) -> tuple[list[int], bool]:
        """v and the positions within minimax distance ``radius`` of it, nearest first, at
        most ``cap`` beyond v; and whether the cap cut it."""
        # ponytail: a Python walk per seed; a Rust kernel (or one walk per tree region) at scale
        out, seen, heap = [], {v}, [(0, v)]
        while heap and len(out) <= cap:
            d, x = heapq.heappop(heap)
            out.append(x)
            lo, hi = self.start[x], self.start[x + 1]
            for u, w in zip(self.dst[lo:hi].tolist(), self.w[lo:hi].tolist(), strict=True):
                if w <= radius and u not in seen:  # a tree: one path to u, so one push
                    seen.add(u)
                    heapq.heappush(heap, (max(d, w), u))
        return out, bool(heap)


def topology_distance(
    identity: NDArray[np.float64] | np.float64, length: NDArray[np.int64], k_t: int
) -> NDArray[np.float64]:
    """Expected |A Δ B| of two proteins' dense k\\_t-mer sets at ``identity``."""
    return 2 * np.maximum(length - k_t + 1, 1) * (1 - identity**k_t)


def held_pairs(
    index: query.Index, hashes: NDArray[np.uint64], seeds: NDArray[np.int64]
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """(seed row, k-mer row) for every k-mer of ``hashes`` each seed position holds."""
    first = np.searchsorted(index.hash, hashes, "left")
    per = np.searchsorted(index.hash, hashes, "right") - first
    rows = np.repeat(first - np.cumsum(per) + per, per) + np.arange(per.sum())
    kmer = np.repeat(np.arange(len(hashes)), per)
    sign = np.asarray(index.sign[rows], dtype=np.int64)
    x = np.concatenate([np.asarray(index.lo[rows]), np.asarray(index.hi[rows])]).astype(np.int64)
    s, kk = np.concatenate([sign, -sign]), np.concatenate([kmer, kmer])
    o = np.lexsort((x, kk))
    x, s, kk = x[o], s[o], kk[o]
    cover = np.cumsum(s)  # each k-mer's signs sum to 0, so the running sum restarts per k-mer
    on = np.flatnonzero((cover[:-1] == 1) & (kk[1:] == kk[:-1]) & (x[1:] > x[:-1]))
    a, b = np.searchsorted(seeds, x[on]), np.searchsorted(seeds, x[on + 1])
    n = b - a
    seed = np.repeat(a - np.cumsum(n) + n, n) + np.arange(n.sum())
    return seed, np.repeat(kk[on], n)


def gather(
    index: query.Index,
    hits: pl.DataFrame,
    seeds: pl.DataFrame,
    pos_len: NDArray[np.uint32],
    min_hits: int,
    tie_z: float,
) -> pl.DataFrame:
    pos = seeds["pos"].to_numpy().astype(np.int64)
    kp, kk = held_pairs(index, hits["hash"].to_numpy(), pos)
    o = np.argsort(kp, kind="stable")
    kp, kk = kp[o], kk[o]
    start = np.searchsorted(kp, np.arange(len(pos) + 1))
    windows = np.maximum(pos_len[pos].astype(np.float64) - index.k + 1, 1)
    scale = windows / np.maximum(seeds["size"].to_numpy().astype(np.float64), 1)
    live = np.ones(hits.height, dtype=bool)
    rows = []
    # ponytail: O(picks x pairs) numpy passes; incremental counts (or _core.gather) at scale
    while True:
        on = live[kk]
        r = np.bincount(kp[on], minlength=len(pos))
        s = r * scale
        best = int(np.argmax(s)) if len(s) else 0
        if not len(s) or r[best] < min_hits:
            break
        mine = np.zeros(len(live), dtype=bool)
        own = kk[start[best] : start[best + 1]]
        mine[own[live[own]]] = True
        shares = np.zeros(len(pos), dtype=bool)
        shares[kp[on & mine[kk]]] = True
        # ponytail: O(nodes) masks per pick for the component and the common subtree
        p = pos[best]
        up = np.flatnonzero((index.node_lo <= p) & (index.node_hi > p))  # best's ancestors
        lo, hi = (index.node_lo[up].min(), index.node_hi[up].max()) if len(up) else (p, p + 1)
        close = s >= s[best] - tie_z * np.sqrt(r[best]) * scale[best]
        tied = np.flatnonzero(shares & close & (pos >= lo) & (pos < hi))  # in best's component
        group = np.isin(kp, tied) & on
        live[kk[group]] = False
        lo, hi, node = int(pos[tied].min()), int(pos[tied].max()) + 1, -1
        if len(tied) > 1:
            span = up[(index.node_lo[up] <= lo) & (index.node_hi[up] >= hi)]
            node = int(span[np.argmin(index.node_hi[span] - index.node_lo[span])])
            lo, hi = int(index.node_lo[node]), int(index.node_hi[node])
        rows.append((len(rows), int(pos[best]), float(s[best]), int(np.unique(kk[group]).size),
                     len(tied), node, lo, hi))  # fmt: skip
    cols = ["rank", "pos", "score", "kmers", "tied", "node", "lo", "hi"]
    out = pl.DataFrame(rows, schema=cols, orient="row")
    return out.join(
        seeds.select("pos", "aai", "aai_lo", "aai_hi"), on="pos", how="left", maintain_order="left"
    )


def ranges(pos: NDArray[np.int64]) -> pl.DataFrame:
    """Sorted distinct positions as half-open ranges."""
    if not len(pos):
        return pl.DataFrame({"lo": pos, "hi": pos})
    cut = np.flatnonzero(np.diff(pos) != 1) + 1
    return pl.DataFrame({"lo": pos[np.r_[0, cut]], "hi": pos[np.r_[cut - 1, len(pos) - 1]] + 1})


def report(path: Path, qdir: Path, a: argparse.Namespace) -> dict:
    index = query.Index(path)
    params = index.meta["params"]
    if params.get("leaves") != "proteins" or float(params.get("topology_t", 0)) != 1:
        raise ValueError("the neighbour rule needs --leaves proteins (dense topology sets)")
    k_t = int(params["topology_k"])
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    leaves = with_aai(pl.read_parquet(qdir / "leaves.parquet"), index.k, a.level)
    nodes = with_aai(pl.read_parquet(qdir / "nodes.parquet"), index.k, a.level)
    leaves.write_parquet(out / "leaves.parquet")
    nodes.write_parquet(out / "nodes.parquet")
    hits = pl.read_parquet(qdir / "hits.parquet")
    pos_len = np.load(path / "pos_len.npy")
    n_pos = len(pos_len)

    # candidates: seeds and their neighbours
    seeds = leaves.filter(pl.col("hit_kmers") >= a.min_hits).sort("pos")
    spos = seeds["pos"].to_numpy().astype(np.int64)
    length = pos_len[spos].astype(np.int64)
    radius = topology_distance(seeds["aai_lo"].to_numpy(), length, k_t) + topology_distance(
        np.float64(1 - a.delta), length, k_t
    )
    tree = Spanning(path, n_pos)
    found, capped = [spos], 0
    for v, r in zip(spos.tolist(), radius.tolist(), strict=True):
        got, cut = tree.cluster(v, r, a.max_neighbours)
        found.append(np.array(got, dtype=np.int64))
        capped += cut
    cand = np.unique(np.concatenate(found))
    candidates = ranges(cand)
    candidates.write_parquet(out / "candidates.parquet")
    members = (
        pl.read_parquet(path / "proteins.parquet")
        .filter(pl.col("pos").is_in(cand))
        .with_columns(seed=pl.col("pos").is_in(spos))
        .sort("pos", "protein_id")
    )
    members.write_parquet(out / "members.parquet")
    candidates_s = time.perf_counter() - t0

    t1 = time.perf_counter()
    gathered = gather(index, hits, seeds, pos_len, a.min_hits, a.tie_z)
    gathered.write_parquet(out / "gather.parquet")
    stats = {
        "seeds": len(spos),
        "candidate_positions": len(cand),
        "neighbour_positions": len(cand) - len(spos),
        "capped_seeds": capped,
        "candidate_ranges": candidates.height,
        "members": members.height,
        "gathered": gathered.height,
        "gathered_tied": int((gathered["tied"] > 1).sum()),
        "gathered_kmers": int(gathered["kmers"].sum()),
        "hit_kmers": hits.height,
        "candidates_s": candidates_s,
        "gather_s": time.perf_counter() - t1,
    }
    (out / "report.json").write_text(json.dumps(stats, indent=1) + "\n")
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--index", required=True)
    parser.add_argument("--query", required=True, help="query.py's --out")
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-hits", type=int, default=1)
    parser.add_argument("--delta", type=float, default=0.1)
    parser.add_argument("--max-neighbours", type=int, default=256)
    parser.add_argument("--tie-z", type=float, default=2.0)
    parser.add_argument("--level", type=float, default=0.95)
    a = parser.parse_args()
    print(json.dumps(report(Path(a.index), Path(a.query), a)))


if __name__ == "__main__":
    main()
