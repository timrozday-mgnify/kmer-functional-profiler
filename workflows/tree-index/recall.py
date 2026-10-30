"""Tree index: candidate recall against spiked truth (docs/tree-index-plan.md, Progress log,
recall-first candidates; the user's target, 2026-10-30: every protein at depth > 0.5 reaches the
candidates, for per-sample tier 2, with few false candidates).

Subcommands:

- ``targets``: ``--n`` proteins with an index position, as the FASTA ``spike.py simulate``
  takes (each protein its own unit).
- ``score``: candidate sets from one query (``query.py``'s ``--out``) under each broadening
  strategy, against ``spike.py simulate``'s ``truth.tsv``. A strain is recalled when its source
  protein's position is a candidate (``unit``: any protein of its MGnify90 cluster is). With
  ``--reads`` and ``--targets``, ``observable`` is the number of the source protein's k-mers
  that occur whole in the strain's own reads (all six frames): an exact k-mer index cannot
  recall the source of a strain with none. Writes
  ``strains.parquet`` (per strain and strategy), ``recall.tsv`` (per strategy, identity and
  depth bin) and ``strategies.tsv`` (per strategy: candidates, and recall at depth > 0.5).

Strategies, each from the leaves with ``min_hits`` hit k-mers (seeds):

- ``leaves``: the seeds only;
- ``gains``: the seeds, and the whole range of every gain posting of a hit k-mer, including
  leaves below that lost it (a gain is posted where a clade first holds a k-mer);
- ``nb_<delta>_<cap>``: report.py's neighbour rule;
- ``component``: every leaf of each seed's tree component (the most any tree walk can add).
"""

import argparse
import gzip
import json
import time
from pathlib import Path

import numpy as np
import polars as pl
import query
import report
from numpy.typing import NDArray
from scipy.sparse import coo_array
from scipy.sparse.csgraph import connected_components

from kmer_functional_profiler import _core
from kmer_functional_profiler.index import read_fasta

DEPTH_BINS = [0.25, 0.5, 1, 2, 5]


def targets(a: argparse.Namespace) -> None:
    members = pl.read_parquet(a.members, columns=["protein_id", "sequence"])
    prot = pl.read_parquet(Path(a.index) / "proteins.parquet").filter(pl.col("pos") >= 0)
    pool = members.join(prot.select("protein_id"), on="protein_id")
    rows = pool.sample(a.n, seed=a.seed).sort("protein_id").rows()
    Path(a.out).write_text("".join(f">{i}\n{s}\n" for i, s in rows))
    print(f"{len(rows)} targets -> {a.out}")


def observable(
    truth: pl.DataFrame, sources: dict[str, str], r1: str, r2: str, k: int
) -> pl.DataFrame:
    """(strain, observable): the source protein's k-mers found in the strain's reads
    (``s<strain>_*``)."""
    reads: dict[int, list[bytes]] = {}
    for path in (r1, r2):
        with gzip.open(path, "rt") as f:
            for name, seq in zip(f, f, strict=False):
                if name.startswith("@s"):
                    reads.setdefault(int(name[2:].split("_")[0]), []).append(seq.strip().encode())
                next(f), next(f)
    rows = []
    for strain, member in truth.select("strain", "member").iter_rows():
        own = _core.hash_proteins([sources[str(member)].encode()], k)["hash"]
        got = reads.get(strain, [])
        seen = _core.hash_dna(got, k, frames="all")["hash"] if got else np.zeros(0, np.uint64)
        rows.append((strain, int(np.isin(own, seen).sum())))
    return pl.DataFrame(rows, schema=["strain", "observable"], orient="row")


def gain_ranges(index: query.Index, hashes: NDArray[np.uint64]) -> NDArray[np.int64]:
    """Every position under a gain posting of ``hashes``."""
    first = np.searchsorted(index.hash, hashes, "left")
    per = np.searchsorted(index.hash, hashes, "right") - first
    rows = np.repeat(first - np.cumsum(per) + per, per) + np.arange(per.sum())
    rows = rows[np.asarray(index.sign[rows]) > 0]
    lo, hi = np.asarray(index.lo[rows], np.int64), np.asarray(index.hi[rows], np.int64)
    n = hi - lo
    return np.unique(np.repeat(lo - np.cumsum(n) + n, n) + np.arange(n.sum()))


def strategies(
    path: Path, qdir: Path, min_hits: int, deltas: list[float], caps: list[int]
) -> dict[str, NDArray[np.int64]]:
    index = query.Index(path)
    k_t = int(index.meta["params"]["topology_k"])
    leaves = report.with_aai(pl.read_parquet(qdir / "leaves.parquet"), index.k, 0.95)
    hits = pl.read_parquet(qdir / "hits.parquet")
    pos_len = np.load(path / "pos_len.npy")
    n_pos = len(pos_len)
    seeds = leaves.filter(pl.col("hit_kmers") >= min_hits).sort("pos")
    spos = seeds["pos"].to_numpy().astype(np.int64)
    out = {"leaves": spos}
    out["gains"] = np.union1d(spos, gain_ranges(index, hits["hash"].to_numpy()))
    tree = report.Spanning(path, n_pos)
    length = pos_len[spos].astype(np.int64)
    d_hit = report.topology_distance(seeds["aai_lo"].to_numpy(), length, k_t)
    for delta in deltas:
        radius = d_hit + report.topology_distance(np.float64(1 - delta), length, k_t)
        for cap in caps:
            got = [tree.cluster(v, r, cap)[0] for v, r in zip(spos, radius.tolist(), strict=True)]
            out[f"nb_{delta}_{cap}"] = np.unique(np.concatenate([spos, *map(np.array, got)]))
    a, b = (np.load(path / f"edge_{x}.npy") for x in "ab")
    graph = coo_array((np.ones(len(a)), (a, b)), shape=(n_pos, n_pos))
    _, comp = connected_components(graph, directed=False)
    out["component"] = np.flatnonzero(np.isin(comp, comp[spos]))
    return out


def score(a: argparse.Namespace) -> None:
    path, out = Path(a.index), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    sets = strategies(path, Path(a.query), a.min_hits, a.deltas, a.caps)
    prot = pl.read_parquet(path / "proteins.parquet")
    truth = (
        pl.read_csv(a.truth, separator="\t", columns=["strain", "member", "identity_bin", "depth"])
        .join(prot.rename({"protein_id": "member"}), on="member", how="left")
        .with_columns(depth_bin=pl.col("depth").cut(DEPTH_BINS, left_closed=True))
    )
    if a.reads:
        sources = dict(read_fasta(a.targets))
        k = int(json.loads((path / "meta.json").read_text())["params"]["k"])
        truth = truth.join(observable(truth, sources, *a.reads, k), on="strain")
    per = truth.with_columns(
        **{s: pl.col("pos").is_in(pl.Series(p).implode()) for s, p in sets.items()},
        **{
            f"unit_{s}": pl.col("cluster_rep").is_in(
                prot.filter(pl.col("pos").is_in(pl.Series(p).implode()))["cluster_rep"].implode()
            )
            for s, p in sets.items()
        },
    )
    per.write_parquet(out / "strains.parquet")
    names = list(sets)
    (
        per.group_by("identity_bin", "depth_bin")
        .agg(n=pl.len(), **{s: pl.col(s).mean() for s in names})
        .sort("identity_bin", "depth_bin", descending=[True, False])
        .with_columns(pl.col("depth_bin").cast(pl.String))
        .write_csv(out / "recall.tsv", separator="\t")
    )
    deep = per.filter(pl.col("depth") > 0.5)
    n_proteins = prot.height
    rows = []
    for s, p in sets.items():
        members = prot.filter(pl.col("pos").is_in(pl.Series(p).implode())).height
        rows.append(
            {
                "strategy": s,
                "positions": len(p),
                "members": members,
                "share_of_index": members / n_proteins,
                "recall_deep": deep[s].mean(),
                "missed_deep": int((~deep[s]).sum()),
                "unit_recall_deep": deep[f"unit_{s}"].mean(),
                "recall_all": per[s].mean(),
            }
        )
    table = pl.DataFrame(rows)
    table.write_csv(out / "strategies.tsv", separator="\t")
    info = {"strains": per.height, "deep": deep.height, "score_s": time.perf_counter() - t0}
    (out / "score.json").write_text(json.dumps(info, indent=1) + "\n")
    with pl.Config(tbl_rows=50, tbl_cols=20, tbl_width_chars=200):
        print(table)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("targets")
    p.add_argument("--members", required=True)
    p.add_argument("--index", required=True)
    p.add_argument("--n", type=int, default=2000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", required=True)
    p = sub.add_parser("score")
    p.add_argument("--index", required=True)
    p.add_argument("--query", required=True)
    p.add_argument("--truth", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--reads", nargs=2, help="the sample's R1 R2, for observable k-mers")
    p.add_argument("--targets", help="targets' FASTA (with --reads)")
    p.add_argument("--min-hits", type=int, default=1)
    p.add_argument("--deltas", type=float, nargs="+", default=[0.0, 0.05, 0.1, 0.15, 0.3])
    p.add_argument("--caps", type=int, nargs="+", default=[16, 256])
    a = parser.parse_args()
    {"targets": targets, "score": score}[a.command](a)


if __name__ == "__main__":
    main()
