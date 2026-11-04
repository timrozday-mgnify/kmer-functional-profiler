"""Tree index, phase 4: benchmark scoring without fixed units (docs/tree-index-plan.md,
Benchmark scoring without fixed units), the same for a tree tier and today's tier 1.

Subcommands:

- ``reported``: what a method reports -> ``--out``: ``units.parquet`` (unit, abundance, aai,
  aai_lo, aai_hi; null when the method has none), ``members.parquet`` (unit, protein_id),
  ``reported.fa``, ``strains.fa`` (the strains' translated CDS) and ``cost.json``. A tree tier
  (``--index``, ``--query``): each candidate leaf at ``--share``, a layer representative's group
  expanded, plus leaves that reads were called to at the read bar; abundance is the reads called
  per residue. Today's tier 1 (``--profile``, from ``kfp query --aai``): each hit unit, expanded
  to its members; abundance is ``coverage_em``.
- ``score``: ``--hits`` (DIAMOND, truth against reported: qseqid, sseqid, pident) -> ``curves.tsv``
  (nearest-match identity curves, recall and precision, overall and by depth and identity bin),
  ``levels.parquet`` (recall and precision down the abundance ranking, per level),
  ``aai.tsv`` (the unit holding each strain's source: estimated AAI against the true identity),
  ``neighbours.tsv`` (with ``--strains``: the neighbour rule's reach on spikes with no own hit)
  and ``score.json`` (read-outs, earth mover's distance, roll-ups).
- ``gate``: pairs of scored runs (same tier name) -> ``gate.tsv``, phase 4's proposed gate.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
import report
from numpy.typing import NDArray
from scipy.optimize import linprog
from scipy.sparse import coo_array

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "spike-benchmark"))
import spike  # noqa: E402

LEVELS = (0.99, 0.95, 0.9, 0.8)
GRID = np.linspace(0.5, 1, 51)
DEPTH_BINS = [0.25, 0.5, 1, 2, 5]
MAX_COST = 1.2  # the gate's query cost and index size, over tier 1's
MAX_SIZE = 1.5


def curve(
    unit: NDArray[np.int64], item: NDArray[np.int64], rank: NDArray[np.int64], n_items: int
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Recall and precision at each cut of a ranking: (unit, item) edges say a reported unit
    matches a true item; ``rank[u]`` is unit u's place in the ranking. Item i is recalled at the
    first unit that matches it."""
    n = len(rank)
    if not n:
        return np.zeros(0), np.zeros(0)
    r = rank[unit]
    first = np.full(n_items, n, dtype=np.int64)
    np.minimum.at(first, item, r)
    found = np.bincount(first[first < n], minlength=n)
    hit = np.zeros(n, dtype=np.int64)
    hit[r] = 1
    return np.cumsum(found) / max(n_items, 1), np.cumsum(hit) / np.arange(1, n + 1)


def at_precision(recall: NDArray, precision: NDArray, p: float) -> float:
    """The most recall at any cut with precision >= p (0 when none reaches it)."""
    ok = precision >= p - 1e-12
    return float(recall[ok].max()) if ok.any() else 0.0


def emd(
    truth: NDArray[np.float64],
    estimate: NDArray[np.float64],
    i: NDArray[np.int64],
    j: NDArray[np.int64],
    identity: NDArray[np.float64],
) -> float:
    """Earth mover's distance between abundances (each normalised to 1), ground distance
    1 - identity, 1 for a pair without a hit. Moving mass along a hit saves its identity, so
    this is 1 - the most identity-weighted flow over hits, a sparse LP.
    ponytail: HiGHS on every hit; a network simplex if a large sample's LP gets slow."""
    if truth.sum() <= 0 or estimate.sum() <= 0 or not len(i):
        return 1.0
    n = np.arange(len(i))
    a = coo_array(
        (np.ones(2 * len(i)), (np.r_[i, len(truth) + j], np.r_[n, n])),
        shape=(len(truth) + len(estimate), len(i)),
    ).tocsr()
    res = linprog(
        -identity,
        A_ub=a,
        b_ub=np.r_[truth / truth.sum(), estimate / estimate.sum()],
        bounds=(0, None),
        method="highs",
    )
    if not res.success:
        raise RuntimeError(f"earth mover's LP: {res.message}")
    return float(1 + res.fun)


def write_fasta(path: Path, rows: Any) -> None:
    with path.open("w") as f:
        for name, seq in rows:
            f.write(f">{name}\n{seq}\n")


def sequences(files: list[str], ids: pl.Series) -> pl.DataFrame:
    return (
        pl.scan_parquet(files)
        .select("protein_id", "sequence", "cluster_rep")
        .filter(pl.col("protein_id").is_in(ids.unique().implode()))
        .collect()
        .unique("protein_id")
    )


def tree_reported(a: argparse.Namespace) -> tuple[pl.DataFrame, pl.DataFrame, dict[str, float]]:
    path, qdir = Path(a.index), Path(a.query)
    prot = pl.read_parquet(path / "proteins.parquet").filter(pl.col("pos") >= 0)
    n_pos = int(prot["pos"].max()) + 1
    gfile = path / "layer_group.npy"
    group = np.load(gfile).astype(np.int64) if gfile.exists() else np.arange(n_pos)
    cand = np.load(qdir / f"cand_verify_{a.share}.npy").astype(np.int64)
    bar = json.loads((qdir / "verify.json").read_text())[f"read_bar_{a.share}"]
    calls = pl.read_parquet(qdir / "reads.parquet").filter(pl.col("score") >= bar)
    if "weight" not in calls.columns:  # tier 1: every tied best leaf counts once
        calls = calls.with_columns(weight=pl.lit(1.0))
    called = calls.group_by("pos").agg(pl.col("weight").sum())
    keep = np.isin(group[:n_pos], cand) | np.isin(np.arange(n_pos), cand)
    keep[called["pos"].to_numpy()] = True
    members = prot.filter(pl.col("pos").is_in(pl.Series(np.flatnonzero(keep)).implode())).select(
        unit="pos", protein_id="protein_id"
    )
    length = (
        pl.scan_parquet(a.members)
        .select("protein_id", length=pl.col("sequence").str.len_chars())
        .join(members.lazy(), on="protein_id")
        .group_by("unit")
        .agg(pl.col("length").max())
        .collect()
    )
    units = (
        length.join(called.rename({"pos": "unit"}), on="unit", how="left")
        .with_columns(abundance=pl.col("weight").fill_null(0) / pl.col("length"))
        .select("unit", "abundance")
    )
    root = qdir.parent if qdir.name == "tier2" else qdir
    aai = root / "leaves.parquet"
    if aai.exists():
        k = int(json.loads((path / "meta.json").read_text())["params"]["k"])
        est = report.with_aai(pl.read_parquet(aai), k, 0.95)
        units = units.join(
            est.select(unit="pos", aai="aai", aai_lo="aai_lo", aai_hi="aai_hi"),
            on="unit",
            how="left",
        )
    t = json.loads((root / "query.json").read_text())["total_s"]
    v = json.loads((root / "verify.json").read_text())
    t += v["sample_total_s"] + v["decoy_total_s"]
    if qdir.name == "tier2":
        t += json.loads((qdir / "tier2.json").read_text())["total_s"]
    return units, members, {"query_s": t}


def flat_reported(a: argparse.Namespace) -> tuple[pl.DataFrame, pl.DataFrame, dict[str, float]]:
    prof = pl.read_csv(a.profile, separator="\t", infer_schema_length=None)
    units = prof.select(
        unit="cluster_rep", abundance=pl.col("coverage_em").cast(pl.Float64).fill_null(0.0),
        aai=pl.col("aai").cast(pl.Float64), aai_lo=pl.col("aai_lo").cast(pl.Float64),
        aai_hi=pl.col("aai_hi").cast(pl.Float64),
    )  # fmt: skip
    members = (
        pl.scan_parquet(a.members)
        .select(unit="cluster_rep", protein_id="protein_id")
        .filter(pl.col("unit").is_in(units["unit"].implode()))
        .unique()
        .collect()
    )
    stages = json.loads(Path(a.stats).read_text())["stages"]
    return units, members, {"query_s": sum(s["wall_s"] for s in stages.values())}


def reported(a: argparse.Namespace) -> None:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    units, members, cost = flat_reported(a) if a.profile else tree_reported(a)
    seqs = sequences(a.members, members["protein_id"])
    units.write_parquet(out / "units.parquet")
    members.write_parquet(out / "members.parquet")
    write_fasta(out / "reported.fa", seqs.select("protein_id", "sequence").iter_rows())
    truth = pl.read_csv(a.truth, separator="\t", columns=["strain", "cds"])
    write_fasta(out / "strains.fa", ((s, spike.translate(c)) for s, c in truth.iter_rows()))
    idx = Path(a.index_dir)
    cost["index_bytes"] = sum(f.stat().st_size for f in idx.rglob("*") if f.is_file())
    (out / "cost.json").write_text(json.dumps(cost, indent=1) + "\n")


def rollup(
    labels: pl.DataFrame, members: pl.DataFrame, truth: pl.DataFrame, units: pl.DataFrame
) -> tuple[NDArray, NDArray]:
    """The curve one level up: a reported label ranks as its most abundant unit, and matches
    when some strain's source holds it."""
    held = members.join(labels, on="protein_id").join(units, on="unit").group_by("label").agg(
        pl.col("abundance").max()
    ).sort("abundance", "label", descending=[True, False])  # fmt: skip
    true = truth.select(protein_id="member").join(labels, on="protein_id")["label"].unique()
    rank = np.arange(held.height, dtype=np.int64)
    item = {x: i for i, x in enumerate(true.to_list())}
    unit_ids, items = [], []
    for u, x in enumerate(held["label"].to_list()):
        if x in item:
            unit_ids.append(u)
            items.append(item[x])
    return curve(
        np.array(unit_ids, dtype=np.int64), np.array(items, dtype=np.int64), rank, len(item)
    )


def shares(x: NDArray[np.float64]) -> list[float]:
    return [float((x >= g - 1e-9).mean()) if len(x) else 0.0 for g in GRID]


def score(a: argparse.Namespace) -> None:
    rep, out = Path(a.reported), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    units = pl.read_parquet(rep / "units.parquet").sort(
        "abundance", "unit", descending=[True, False]
    )
    members = pl.read_parquet(rep / "members.parquet")
    truth = (
        pl.read_csv(a.truth, separator="\t", columns=["strain", "member", "identity_bin",
                                                      "identity", "depth"])
        .sort("strain")
        .with_columns(depth_bin=pl.col("depth").cut(DEPTH_BINS, left_closed=True).cast(pl.String))
    )  # fmt: skip
    hits = pl.read_csv(
        a.hits, separator="\t", has_header=False, new_columns=["strain", "protein_id", "pident"]
    )
    rank = units.with_row_index("rank").select("unit", pl.col("rank").cast(pl.Int64))
    n_strains = truth.height
    strain_i = truth.with_row_index("i").select("strain", pl.col("i").cast(pl.Int64))
    edges = (
        hits.join(members, on="protein_id")
        .group_by("strain", "unit")
        .agg(identity=pl.col("pident").max() / 100)
        .join(rank, on="unit")
        .join(strain_i, on="strain")
    )
    i = edges["i"].to_numpy()
    r = edges["rank"].to_numpy()
    identity = edges["identity"].to_numpy()

    # nearest-match curves: each strain's nearest unit; each unit's nearest strain
    nearest_t = np.zeros(n_strains)
    np.maximum.at(nearest_t, i, identity)
    nearest_u = np.zeros(units.height)
    np.maximum.at(nearest_u, r, identity)
    rows = [("recall", "all", "all", nearest_t), ("precision", "all", "all", nearest_u)]
    for by in ("identity_bin", "depth_bin"):
        col = truth[by].cast(pl.String).to_numpy()
        rows += [("recall", by, v, nearest_t[col == v]) for v in np.unique(col.astype(str))]
    pl.DataFrame(
        [
            {
                "side": s,
                "by": by,
                "bin": v,
                "n": len(x),
                **dict(zip(map(str, GRID.round(2)), shares(x), strict=True)),
            }
            for s, by, v, x in rows
        ]  # fmt: skip
    ).write_csv(out / "curves.tsv", separator="\t")

    # read-outs, and each level's recall and precision down the ranking
    info: dict[str, Any] = {"strains": n_strains, "units": units.height}
    levels = []
    for lev in LEVELS:
        ok = identity >= lev - 1e-9
        rec, prec = curve(r[ok], i[ok], np.arange(units.height, dtype=np.int64), n_strains)
        info[f"recall_{lev}"] = float(rec[-1]) if len(rec) else 0.0
        info[f"precision_{lev}"] = float(prec[-1]) if len(prec) else 0.0
        levels.append(pl.DataFrame({"level": str(lev), "recall": rec, "precision": prec}))

    # abundance; a unit without reads holds none
    ab = units["abundance"].fill_null(0).to_numpy()
    info["emd"] = emd(truth["depth"].to_numpy(), ab, i, r, identity)

    # roll-ups: MGnify90 through membership, and any --rollup name=table
    tables = {
        "mgnify90": pl.scan_parquet(a.members)
        .select("protein_id", label="cluster_rep")
        .filter(
            pl.col("protein_id").is_in(
                pl.concat([members["protein_id"], truth["member"]]).unique().implode()
            )
        )
        .collect()
    }
    for spec in a.rollup:
        name, file = spec.split("=", 1)
        t = pl.read_parquet(file)
        tables[name] = t.select(pl.col(t.columns[0]).alias("protein_id"),
                                label=pl.col(t.columns[1]).cast(pl.String)).unique()  # fmt: skip
    for name, labels in tables.items():
        rec, prec = rollup(labels, members, truth, units)
        info[f"recall_{name}"] = float(rec[-1]) if len(rec) else 0.0
        info[f"precision_{name}"] = float(prec[-1]) if len(prec) else 0.0
        levels.append(pl.DataFrame({"level": name, "recall": rec, "precision": prec}))
    pl.concat(levels).write_parquet(out / "levels.parquet")

    # AAI calibration: the unit holding each strain's source
    if "aai" in units.columns:
        cal = (
            truth.join(members.rename({"protein_id": "member"}), on="member")
            .join(units, on="unit")
            .filter(pl.col("aai").is_not_null())
            .group_by("strain")
            .first()
            .group_by("identity_bin")
            .agg(
                n=pl.len(),
                aai=pl.col("aai").median(),
                bias=(pl.col("aai") - pl.col("identity")).median(),
                covered=pl.col("identity").is_between(pl.col("aai_lo"), pl.col("aai_hi")).mean(),
            )
            .sort("identity_bin", descending=True)
        )
        cal.write_csv(out / "aai.tsv", separator="\t")

    # the neighbour rule: spikes without an own hit that it still reaches (recall.py score)
    if a.strains:
        per = pl.read_parquet(a.strains).filter(~pl.col("leaves"))
        rules = [c for c in per.columns if c.startswith("nb_")]
        info["unhit_spikes"] = per.height
        per.with_columns(pl.col("depth_bin").cast(pl.String)).group_by(
            "identity_bin", "depth_bin"
        ).agg(n=pl.len(), **{c: pl.col(c).mean() for c in rules}).sort(
            "identity_bin", "depth_bin", descending=[True, False]
        ).write_csv(out / "neighbours.tsv", separator="\t")
    info |= json.loads((rep / "cost.json").read_text())
    (out / "score.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info, indent=1))


def gate(a: argparse.Namespace) -> None:
    runs = {Path(d).name: Path(d) for d in a.runs}
    rows = []
    for name, tree in sorted(runs.items()):
        if not name.startswith("tree-"):
            continue
        flat = runs.get("flat-" + name.removeprefix("tree-").removesuffix("-tier2"))
        if flat is None:
            continue
        ti, fi = (json.loads((d / "score.json").read_text()) for d in (tree, flat))
        tl, fl = (pl.read_parquet(d / "levels.parquet") for d in (tree, flat))
        for lev in fl["level"].unique(maintain_order=True).to_list():
            t = tl.filter(pl.col("level") == lev)
            fp, fr = fi[f"precision_{lev}"], fi[f"recall_{lev}"]
            matched = at_precision(t["recall"].to_numpy(), t["precision"].to_numpy(), fp)
            rows.append(
                {"tree": name, "measure": lev, "tier1": fr, "tree_value": matched,
                 "tier1_precision": fp, "bar": fr, "passes": matched >= fr}
            )  # fmt: skip
        rows += [
            {"tree": name, "measure": "query_s", "tier1": fi["query_s"],
             "tree_value": ti["query_s"], "tier1_precision": None,
             "bar": MAX_COST * fi["query_s"], "passes": ti["query_s"] <= MAX_COST * fi["query_s"]},
            {"tree": name, "measure": "index_bytes", "tier1": fi["index_bytes"],
             "tree_value": ti["index_bytes"], "tier1_precision": None,
             "bar": MAX_SIZE * fi["index_bytes"],
             "passes": ti["index_bytes"] <= MAX_SIZE * fi["index_bytes"]},
            {"tree": name, "measure": "emd (reported, not gated)", "tier1": fi["emd"],
             "tree_value": ti["emd"], "tier1_precision": None, "bar": None, "passes": None},
        ]  # fmt: skip
    table = pl.DataFrame(rows, infer_schema_length=None)
    table.write_csv(a.out, separator="\t")
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=160):
        print(table)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("reported")
    p.add_argument("--truth", required=True)
    p.add_argument("--members", nargs="+", required=True)
    p.add_argument("--index-dir", required=True, help="the index measured for size")
    p.add_argument("--out", required=True)
    p.add_argument("--index", help="tree tier: its index")
    p.add_argument("--query", help="tree tier: query.py's output, or its tier2/")
    p.add_argument("--share", type=float, default=0.01)
    p.add_argument("--profile", help="today's tier 1: kfp query --aai's TSV")
    p.add_argument("--stats", help="today's tier 1: kfp query --stats")
    p.set_defaults(func=reported)
    p = sub.add_parser("score")
    p.add_argument("--reported", required=True)
    p.add_argument("--truth", required=True)
    p.add_argument("--hits", required=True)
    p.add_argument("--members", nargs="+", required=True)
    p.add_argument("--rollup", nargs="*", default=[], help="name=parquet (protein_id, label)")
    p.add_argument("--strains", help="recall.py score's strains.parquet")
    p.add_argument("--out", required=True)
    p.set_defaults(func=score)
    p = sub.add_parser("gate")
    p.add_argument("runs", nargs="+", help="score dirs: tree-<tier>[-tier2], flat-<tier>")
    p.add_argument("--out", required=True)
    p.set_defaults(func=gate)
    a = parser.parse_args()
    a.func(a)


if __name__ == "__main__":
    sys.exit(main())
