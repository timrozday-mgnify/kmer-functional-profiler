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
- ``reads``: precision and read classification from ``verify.py``'s output. Candidates
  (``cand_verify_<share>``, groups expanded): the share that are a source protein's leaf or
  in a source's MGnify90 unit (``cluster_rep``), and of their units, the share that are a
  source's. Reads: each read pair's best leaves (``reads.parquet``, ties shared equally, or by
  its ``weight`` column) at the
  read-level bar (``verify.json``'s ``read_bar_<share>``), against its strain (read names
  ``s<strain>_<n>``; others, ``spike.py``'s background, counted apart): sensitivity (called
  to the source's unit, over the strain's reads) and precision by leaf, layer group, unit
  and, with ``--pfam``, Pfam family (the called leaf's and the source's accessions overlap).
  Writes ``reads.tsv`` (per identity bin) and ``reads.json``.
- ``summary``: one row per build: postings, decoy bar and candidate recall by identity; tier 1
  as a profile (its candidate units' precision and recall, strains detected at depth >= 1)
  and, from ``tier2/``, tier 2 as a read classifier (read recall and precision, leaf, top-leaf
  and Pfam precision, median hits, its units' precision and recall), across builds.

Strategies, each from the leaves with ``min_hits`` hit k-mers (seeds):

- ``leaves``: the seeds only;
- ``gains``: the seeds, and the whole range of every gain posting of a hit k-mer, including
  leaves below that lost it (a gain is posted where a clade first holds a k-mer);
- ``nb_<delta>_<cap>``: report.py's neighbour rule;
- ``component``: every leaf of each seed's tree component (the most any tree walk can add);
- ``ev_<share>`` (and ``_nb``, with the last delta and cap's neighbour rule): seeds by
  evidence from one read pair (query.py's ``seeds.parquet``), at the least bits whose expected
  chance seeds are at most ``share`` of the seeds.
"""

import argparse
import gzip
import json
import time
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
import query
import report
from numpy.typing import NDArray
from polars.testing import assert_frame_equal
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
    path: Path,
    qdir: Path,
    min_hits: int,
    deltas: list[float],
    caps: list[int],
    shares: tuple[float, ...] = (0.001, 0.01, 0.1),
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
    aai_lo = dict(zip(leaves["pos"].to_list(), leaves["aai_lo"].to_list(), strict=True))

    def walk(start: NDArray[np.int64], delta: float, cap: int) -> NDArray[np.int64]:
        length = pos_len[start].astype(np.int64)
        lo = np.array([aai_lo.get(int(v), 0.0) for v in start])
        radius = report.topology_distance(lo, length, k_t) + report.topology_distance(
            np.float64(1 - delta), length, k_t
        )
        got = [tree.cluster(v, r, cap)[0] for v, r in zip(start, radius.tolist(), strict=True)]
        return np.unique(np.concatenate([start, *map(np.array, got)]))

    for delta in deltas:
        for cap in caps:
            out[f"nb_{delta}_{cap}"] = walk(spos, delta, cap)
    # seeds by evidence (query.py --reads with key bits), at the bar for each chance share
    if (qdir / "seeds.parquet").exists():
        ev = pl.read_parquet(qdir / "seeds.parquet")
        curve = pl.read_parquet(qdir / "evidence.parquet")
        for share in shares:
            ok = curve.filter(pl.col("chance") <= share * pl.col("seeds").clip(lower_bound=1))
            bar = ok["bits"].min() if ok.height else curve["bits"].max()
            start = ev.filter(pl.col("score") >= bar)["pos"].to_numpy().astype(np.int64)
            out[f"ev_{share}"] = start
            out[f"ev_{share}_nb"] = walk(start, deltas[-1], caps[-1])
    for f in sorted(qdir.glob("cand_*.npy")):  # query.py --decoy's candidates
        start = np.load(f).astype(np.int64)
        out[f.stem[5:]] = start
        out[f"{f.stem[5:]}_nb"] = walk(start, deltas[-1], caps[-1])
    if (path / "layer_group.npy").exists():  # a layer on representatives: each seeds its group
        group = np.load(path / "layer_group.npy").astype(np.int64)
        o = np.argsort(group, kind="stable")
        first = np.searchsorted(group[o], np.arange(len(group) + 1))
        reps = lambda p: p[group[p] == p]  # noqa: E731
        out = {s: np.union1d(p, np.concatenate([o[first[r] : first[r + 1]] for r in reps(p)]
                                                 or [p[:0]])) for s, p in out.items()}  # fmt: skip
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


MGNIFY_PROTEINS = 5_739_566_171  # MGnify Proteins 2026_07 (docs/mgnify-2026_07.md)
MGNIFY_MEAN_LENGTH = 1.07e12 / MGNIFY_PROTEINS  # residues per protein, ~186


def release_gb(idx: Path, meta: dict, n_release: int = MGNIFY_PROTEINS) -> dict[str, float]:
    """Tier 1's footprint scaled to ``n_release`` proteins by this index's per-protein rates
    (docs/tree-index-plan.md, Progress log: size estimates revisited), in two parts:

    - *lookup*: postings (hash 8 B, lo and hi 4 B each, sign 1 B, bits 1 B, k code 1-2 B),
      ``pos_size`` and the node tables (lo, hi, level, size), and the protein -> leaf -> unit
      map (``proteins.parquet`` as written);
    - *align*: what ``verify.py`` aligns against, to be reconsidered with the rest of the
      index: ``layer_group`` and each group's representative (or every leaf without groups;
      residues 1 B at the release's mean length, offsets 8 B).

    Positions stay u32 with the release in two shards. Rates follow density: a slice's
    (complete neighbourhoods) is the release's best guess, a cluster sample's an upper bound.
    Tier 2's full store and the spanning-tree edges are apart. ponytail: postings are not
    rescaled to the release's mean length; a long-protein sample overstates the layer's."""
    s, f = meta["stats"], meta["files"]
    n = s["proteins"]
    lengths = np.load(idx / "pos_len.npy").astype(np.int64)
    groups = idx / "layer_group.npy"
    reps = np.unique(np.load(groups)) if groups.exists() else np.arange(len(lengths))
    nodes = (f["node_lo.npy"] - 128) // 4  # u32 after the .npy header
    code = (f["kcode.npy"] - 128) / max(s["postings"], 1) if "kcode.npy" in f else 1
    residues = lengths[reps].sum() * MGNIFY_MEAN_LENGTH / lengths.mean()
    per = {
        "postings": (18 + code) * s["postings"] / n,
        "tables": (4 * s["positions"] + 16 * nodes + f["proteins.parquet"]) / n,
        "align": (4 * groups.exists() * s["positions"] + residues + 8 * len(reps)) / n,
    }
    got = {f"release_{k}_gb": round(v * n_release / 1e9) for k, v in per.items()}
    lookup = got["release_postings_gb"] + got["release_tables_gb"]
    return got | {
        "release_lookup_gb": lookup,
        "release_tier1_gb": lookup + got["release_align_gb"],
    }


def release_pairs(q: dict, verify: dict, n: int, n_release: int = MGNIFY_PROTEINS) -> dict:
    """``verify.py``'s alignment work per read pair, now and scaled to ``n_release`` proteins:
    chance pairs (the decoy reads') grow with the database's keys, so with its proteins; true
    pairs do not. ponytail: linear in proteins, as chance seeds are; a denser release shares
    more keys, so this is an upper bound."""
    reads = max(q["read_pairs"], 1)
    sample, decoy = verify["sample_pairs"] / reads, verify["decoy_pairs"] / reads
    return {
        "pairs_per_read": round(sample),
        "chance_share": round(decoy / max(sample, 1e-9), 2),
        "release_pairs_per_read": round(max(sample - decoy, 0) + decoy * n_release / n),
    }


def summary(a: argparse.Namespace) -> None:
    """One row per build under ``--root`` (``idx-<name>`` or ``index-<name>``, ``q-<name>``,
    ``s-<name>``): its
    layers, postings (and over ``--base``'s), the query's bar and decoy seeds, and recall at
    depth > ``--depth`` per identity for ``--strategy``, with its candidates."""
    root = Path(a.root)
    rows = []
    for sdir in sorted(root.glob("s-*")):
        name = sdir.name[2:]
        idx = root / f"idx-{name}"
        idx = idx if idx.exists() else root / f"index-{name}"
        meta = json.loads((idx / "meta.json").read_text())
        q = json.loads((root / f"q-{name}" / "query.json").read_text())
        per = pl.read_parquet(sdir / "strains.parquet").filter(pl.col("depth") > a.depth)
        cand = pl.read_csv(sdir / "strategies.tsv", separator="\t")
        row = {
            "name": name,
            "layers": meta["params"].get("layers") or "",
            "postings": meta["stats"]["postings"],
            "bytes": meta["stats"]["bytes"],
            "bar": q.get("evidence_bits"),
            "decoy_seeds": q.get("decoy_seeds"),
            "query_s": round(q["total_s"]),
            "candidates": cand.filter(pl.col("strategy") == a.strategy)["positions"].item(),
        } | release_gb(idx, meta)
        vfile = root / f"q-{name}" / "verify.json"
        if vfile.exists():
            row |= release_pairs(q, json.loads(vfile.read_text()), meta["stats"]["proteins"])
        # tier 1 is a profile: its candidates' units and strain detection; tier 2 classifies
        # reads (tier 1's read numbers stay in its reads.json, a diagnostic)
        if (sdir / "reads.json").exists():
            r = json.loads((sdir / "reads.json").read_text())
            row |= {k: r.get(k) for k in ("cand_unit_precision", "unit_recall", "detected_depth1")}
        if (sdir / "tier2" / "reads.json").exists():
            r = json.loads((sdir / "tier2" / "reads.json").read_text())
            keys = (
                "read_recall", "read_precision", "leaf_precision", "leaf_precision_top",
                "pfam_precision", "hits_median", "cand_unit_precision", "unit_recall",
            )  # fmt: skip
            row |= {f"t2_{k}": r.get(k) for k in keys}
        for ident, got in per.group_by("identity_bin").agg(pl.col(a.strategy).mean()).iter_rows():
            row[f"r{ident}"] = got
        rows.append(row)
    table = pl.DataFrame(rows)
    if a.base in table["name"]:
        base = table.filter(pl.col("name") == a.base)["postings"].item()
        table = table.with_columns(x_base=pl.col("postings") / base)
    recall = sorted((c for c in table.columns if c[0] == "r" and c[1:2].isdigit()), reverse=True)
    table = table.select(*(c for c in table.columns if c not in recall), *recall).sort(recall[-1])
    table.write_csv(Path(a.out), separator="\t", float_precision=3)
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=250, float_precision=3):
        print(table)


DEPTH_EDGES = [0.5, 1, 2, 5]
CODING_EDGES = [0.25, 0.5, 0.75]  # the read pair's share of bases in the coding sequence
MIN_STRAINS = 30  # a cell of strain means: bins merged below this


def mean_ci(x: NDArray[np.float64], seed: int = 0, n: int = 1000) -> tuple[float, float, float]:
    """Mean and a 95% bootstrap interval (resampling the values: strains)."""
    if not len(x):
        return float("nan"), float("nan"), float("nan")
    boot = x[np.random.default_rng(seed).integers(0, len(x), (n, len(x)))].mean(axis=1)
    return float(x.mean()), float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def merge_bins(counts: list[int], least: int) -> list[int]:
    """Each bin's merged bin: adjacent bins pooled in order until each pool holds ``least``
    (a short tail joins the pool before it)."""
    out, pool, held = [], 0, 0
    for c in counts:
        out.append(pool)
        held += c
        if held >= least:
            pool, held = pool + 1, 0
    if held and pool:  # the tail is short: into the last full pool
        out = [min(o, pool - 1) for o in out]
    return out


def span(pools: list[int], pool: int) -> str:
    """A pool of depth bins as a range, from DEPTH_EDGES (``0-0.5``, ``5-``)."""
    edges = [0, *DEPTH_EDGES, None]
    mine = [i for i, q in enumerate(pools) if q == pool]
    hi = edges[mine[-1] + 1]
    return f"{edges[mine[0]]}-{'' if hi is None else hi}"


def reads(a: argparse.Namespace) -> None:
    """Read classification and presence, scored so each number reads directly (Progress log:
    an intuitive scoring). A read pair is *called* when its best alignment clears the
    read-level bar; it is *right* (unit, leaf, Pfam) when any of its tied best leaves is the
    source's. Read recall is right calls over the strains' read pairs, read precision right
    calls over the strains' called pairs (the background's calls apart: their truth is
    unknown). Per-strain means with 95% bootstrap intervals carry the variation between
    strains (cells of strain means pool depth bins below MIN_STRAINS strains); strain detection
    is a strain with a right call."""
    path, qdir, out = Path(a.index), Path(a.query), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    prot = pl.read_parquet(path / "proteins.parquet").filter(pl.col("pos") >= 0)
    n_pos = int(prot["pos"].max()) + 1
    gfile = path / "layer_group.npy"
    group = np.load(gfile).astype(np.int64) if gfile.exists() else np.arange(n_pos)
    truth = pl.read_csv(
        a.truth, separator="\t", columns=["strain", "member", "identity_bin", "depth"]
    ).join(
        prot.select(member="protein_id", src_pos="pos", unit="cluster_rep"),
        on="member",
        how="left",
    )
    units = prot.group_by("pos").agg(units=pl.col("cluster_rep").unique())
    info: dict[str, Any] = {}
    # candidates, each representative's group expanded
    cand = np.load(qdir / f"cand_verify_{a.share}.npy").astype(np.int64)
    cand = np.flatnonzero(np.isin(group[:n_pos], cand) | np.isin(np.arange(n_pos), cand))
    cu = units.filter(pl.col("pos").is_in(pl.Series(cand).implode()))
    true_units = set(truth["unit"].drop_nulls().to_list())
    cand_units = set(cu.explode("units")["units"].to_list())
    info |= {
        "candidates": len(cand),
        "cand_leaf_precision": float(
            np.isin(cand, truth["src_pos"].drop_nulls().to_numpy()).mean()
        ),
        "cand_units": len(cand_units),
        "cand_unit_precision": len(cand_units & true_units) / max(len(cand_units), 1),
        "unit_recall": len(cand_units & true_units) / max(len(true_units), 1),
    }
    # read pairs in file order: the simulator's reads.parquet, else the names (s<strain>_<n>)
    if a.read_truth:
        pairs = pl.read_parquet(a.read_truth)
        pairs = pairs.with_columns(coding=pl.col("coding") / pl.col("pair_bases"))
    else:
        with gzip.open(a.reads, "rt") as f:
            names = [
                int(line[2:].split("_")[0]) if line.startswith("@s") else -1
                for i, line in enumerate(f) if i % 4 == 0
            ]  # fmt: skip
        pairs = pl.DataFrame({"read": np.arange(len(names)), "strain": names, "coding": None})
    pairs = pairs.select("read", "strain", pl.col("coding").cast(pl.Float64))
    bar = json.loads((qdir / "verify.json").read_text())[f"read_bar_{a.share}"]
    calls = pl.read_parquet(qdir / "reads.parquet").filter(pl.col("score") >= bar)
    calls = calls.join(pairs.select("read", "strain"), on="read")
    bg_called = calls.filter(pl.col("strain") < 0)["read"].n_unique()
    calls = calls.join(truth, on="strain").join(units, on="pos", how="left")
    if "weight" not in calls.columns:  # tier 1: no shares, its first tie is its top
        calls = calls.with_columns(weight=pl.lit(1.0))
    sp = calls["src_pos"].fill_null(-1).to_numpy()
    # per hit, right at each level (``*_ok``: apart from the truth's own columns)
    calls = calls.with_columns(
        leaf_ok=(pl.col("pos") == pl.col("src_pos")).fill_null(False),
        group_ok=pl.Series(
            (sp >= 0) & (group[calls["pos"].to_numpy()] == group[np.maximum(sp, 0)])
        ),
        unit_ok=pl.col("units").list.contains(pl.col("unit")).fill_null(False),
    )
    levels = ["unit_ok", "group_ok", "leaf_ok"]
    if a.pfam:
        fam = pl.read_parquet(a.pfam, columns=["protein_id", "pfam_accession"]).unique()
        leaf_fam = prot.join(fam, on="protein_id").group_by("pos").agg(lf=pl.col("pfam_accession"))
        src_fam = fam.group_by(member="protein_id").agg(sf=pl.col("pfam_accession"))
        calls = calls.join(leaf_fam, on="pos", how="left").join(src_fam, on="member", how="left")
        calls = calls.with_columns(
            pfam_ok=pl.when(pl.col("sf").is_null())
            .then(None)
            .otherwise(pl.col("lf").list.set_intersection("sf").list.len().fill_null(0) > 0)
        )
        levels.append("pfam_ok")
        # a source without a Pfam hit has no family to agree with: its reads are left out of
        # Pfam precision (before, its null turned into a wrong call)
        annotated = src_fam["member"].implode()
    per = (
        calls.sort("weight", descending=True)
        .group_by("read")
        .agg(
            hits=pl.len(),
            hit_units=pl.col("units").flatten().n_unique(),
            leaf_top=pl.col("leaf_ok").first(),
            **{x: pl.col(x).any() for x in levels},
        )
    )
    allr = (
        pairs.filter(pl.col("strain") >= 0)
        .join(truth, on="strain")
        .join(per, on="read", how="left")
        .with_columns(called=pl.col("hits").is_not_null(), **{x: pl.col(x).fill_null(False)
                                                              for x in (*levels, "leaf_top")})
    )  # fmt: skip
    if a.pfam:
        allr = allr.with_columns(pfam_known=pl.col("member").is_in(annotated))
    called = allr.filter(pl.col("called"))
    # every strain, a strain with no read pair undetected (its recall has no reads: null)
    strains = (
        truth.select("strain", "identity_bin", "depth")
        .join(
            allr.group_by("strain").agg(
                recall=pl.col("unit_ok").mean(), detected=pl.col("unit_ok").any()
            ),
            on="strain",
            how="left",
        )
        .with_columns(pl.col("detected").fill_null(False))
    )

    def rates(df: pl.DataFrame) -> dict[str, Any]:
        c = df.filter(pl.col("called"))
        got = {
            "reads": df.height,
            "called": c.height,
            "read_recall": df["unit_ok"].mean(),
            "read_precision": c["unit_ok"].mean(),
            "group_precision": c["group_ok"].mean(),
            "leaf_recall": df["leaf_ok"].mean(),
            "leaf_precision": c["leaf_ok"].mean(),
            "leaf_precision_top": c["leaf_top"].mean(),
            "hits_median": c["hits"].median(),
            "units_median": c["hit_units"].median(),
        }
        if "pfam_ok" in levels:
            known = c.filter(pl.col("pfam_known"))
            got["pfam_precision"] = known["pfam_ok"].mean()
            got["pfam_known_share"] = known.height / c.height if c.height else None
        return got

    rows = []
    for ident, df in sorted(allr.partition_by("identity_bin", as_dict=True).items(), reverse=True):
        m, lo, hi = mean_ci(
            strains.filter(pl.col("identity_bin") == ident[0])["recall"].drop_nulls().to_numpy()
        )  # fmt: skip
        rows.append({"identity_bin": ident[0], **rates(df), "strains": strains.filter(
            (pl.col("identity_bin") == ident[0]) & pl.col("recall").is_not_null()).height,
            "read_recall_strain": m,
            "ci_lo": lo, "ci_hi": hi})  # fmt: skip
    pl.DataFrame(rows).write_csv(out / "reads.tsv", separator="\t", float_precision=4)
    # strain detection and strain-mean recall by depth x identity, depth bins pooled per row
    labels = ["<0.5", "0.5-1", "1-2", "2-5", ">=5"]
    strains = strains.with_columns(
        d=pl.col("depth").cut(DEPTH_EDGES, labels=labels, left_closed=True)
    )
    det = []
    for ident in sorted(strains["identity_bin"].unique().to_list(), reverse=True):
        s_i = strains.filter(pl.col("identity_bin") == ident)
        counts = [s_i.filter(pl.col("d") == b).height for b in labels]
        pools = merge_bins(counts, MIN_STRAINS)
        for p_ in sorted(set(pools)):
            names = [b for b, q in zip(labels, pools, strict=True) if q == p_]
            cell = s_i.filter(pl.col("d").cast(pl.String).is_in(names))
            d_m, d_lo, d_hi = mean_ci(cell["detected"].cast(pl.Float64).to_numpy())
            r_m, r_lo, r_hi = mean_ci(cell["recall"].drop_nulls().to_numpy())
            det.append({
                "identity_bin": ident, "depth": span(pools, p_), "strains": cell.height,
                "detected": d_m, "detected_lo": d_lo, "detected_hi": d_hi,
                "read_recall_strain": r_m, "recall_lo": r_lo, "recall_hi": r_hi,
            })  # fmt: skip
    pl.DataFrame(det).write_csv(out / "detection.tsv", separator="\t", float_precision=3)
    if a.read_truth:  # by the pair's coding share: what a read can show
        shares = ["<25%", "25-50%", "50-75%", ">=75%"]
        by = allr.with_columns(share=pl.col("coding").cut(CODING_EDGES, labels=shares))
        parts = by.partition_by("share", as_dict=True)
        rows = [{"coding": x, **rates(parts[(x,)])} for x in shares if (x,) in parts]
        pl.DataFrame(rows).write_csv(out / "reads_by_coding.tsv", separator="\t", float_precision=4)
    deep = strains.filter(pl.col("depth") >= 1)
    m, lo, hi = mean_ci(strains["recall"].drop_nulls().to_numpy())
    info |= {
        "read_bar": bar,
        **rates(allr),
        "read_precision_bg_lower": called["unit_ok"].sum() / max(called.height + bg_called, 1),
        "read_recall_strain": m,
        "read_recall_strain_ci": [lo, hi],
        "detected_depth1": float(deep["detected"].mean()) if deep.height else None,
        "strains": strains.height,
        "background_reads": int((pairs["strain"] < 0).sum()),
        "background_called": bg_called,
    }
    (out / "reads.json").write_text(json.dumps(info, indent=1, default=float) + "\n")
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=220, float_precision=3):
        print(pl.read_csv(out / "reads.tsv", separator="\t"))
        print(pl.read_csv(out / "detection.tsv", separator="\t"))
    print(json.dumps(info, default=float))


def chunks(a: argparse.Namespace) -> None:
    """Whether tier 2's runs at several memory budgets (``tier2.py --memory-gb``) made the same
    calls and hand-off assignments: exactly, and to float rounding (EM weights)."""
    runs = sorted(Path(r) for r in a.runs)
    files = ("reads.parquet", "handoff/assignments.parquet")
    read = lambda run, f: (lambda d: d.sort(d.columns))(pl.read_parquet(run / f))  # noqa: E731
    same: dict[str, dict[str, bool]] = {}
    for f in files:
        first = read(runs[0], f)
        exact, close = True, True
        for run in runs[1:]:
            other = read(run, f)
            exact &= other.equals(first)
            try:
                assert_frame_equal(other, first, check_exact=False, rtol=1e-5)
            except AssertionError:
                close = False
        same[f] = {"exact": exact, "close": close}
    info = {r.name: json.loads((r / "tier2.json").read_text()) for r in runs}
    Path(a.out).write_text(json.dumps({"identical": same, "runs": info}, indent=1) + "\n")
    print(json.dumps(same))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("targets")
    p.add_argument("--members", nargs="+", required=True)
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
    p = sub.add_parser("reads")
    p.add_argument("--index", required=True)
    p.add_argument("--query", required=True, help="query.py and verify.py's --out")
    p.add_argument("--truth", required=True)
    p.add_argument("--reads", help="the sample's R1 (read names; without --read-truth)")
    p.add_argument("--read-truth", help="spike.py's reads.parquet: strain, coding per pair")
    p.add_argument("--pfam", help="protein_id, pfam_accession parquet (Pfam precision)")
    p.add_argument("--share", default="0.01", help="the verify.py bar's share")
    p.add_argument("--out", required=True)
    p = sub.add_parser("summary")
    p.add_argument("--root", required=True, help="holds idx-<name>, q-<name> and s-<name>")
    p.add_argument("--strategy", default="ev_0.01_nb")
    p.add_argument("--depth", type=float, default=1.0)
    p.add_argument("--base", default="S3", help="postings are also given over this build's")
    p.add_argument("--out", required=True)
    p = sub.add_parser("chunks")
    p.add_argument("--runs", nargs="+", required=True, help="tier2.py --out directories")
    p.add_argument("--out", required=True)
    a = parser.parse_args()
    {"targets": targets, "score": score, "reads": reads, "summary": summary, "chunks": chunks}[
        a.command
    ](a)


if __name__ == "__main__":
    main()
