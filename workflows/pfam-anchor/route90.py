"""The tree index's route on its own, with no anchoring, on the compact index
(``workflows/tree-index/compact.py``): tier 1 on MGnify90 units (~90% AAI), each unit's keys on
its representative (MGnify's cluster_rep), k 9, about 1 in 32 of the representatives' residues
keyed and spent per unit in proportion to sqrt(members), a read-pair floor of 30 bits; each
candidate (read pair, unit) pseudo-aligned to the representative to confirm the unit; then every
read pair routed to the confirmed representatives and pseudo-aligned (no leaves, no full
alignment). Recall and precision by depth x AAI (docs/tree-index-plan.md, Progress log: the route
on its own, 2026-11-09; the compact route, 2026-11-10).

Steps, into ``--out``:

- ``reads``: ``strains/`` (reference proteins mutated to 100-70% identity, BLOSUM62) and
  ``relatives/`` (pfam-anchor's held-out proteins, at their real identity to the nearest
  protein left in the reference), each simulated at every depth of :data:`DEPTHS` (one
  ``spike.py simulate`` per depth, at that depth exactly) and merged; ``truth.parquet`` per set:
  strain, depth, aai, best (the indexed protein a perfect method would report).
- ``build``: the compact index on pfam-anchor's reference (``compact.py shard`` per
  ``--shards``, then ``merge``) into ``index-<tag>``, the tag naming its sparsity (``r32``, the
  standard index, or ``r8`` and ``r4``, the larger ones: :data:`RATES`).
- ``query``: per set, ``compact.py query`` into ``<set>/c-<tag>`` (tier 1 and confirmation: the
  index's cost), then ``recruit.py`` into ``<set>/r-<tag>`` (recruitment to the confirmed
  representatives, a step of its own, timed apart); bars at ``--chance-scale`` (the release's
  proteins over the reference's, by default).
- ``score``: ``route90_recall_<tag>.tsv`` (per set x depth x AAI: target recall with Wilson
  intervals over strains, within 0.7 and 0.5, the unit's, read recall and precision) and
  ``route90_targets_<tag>.tsv`` (present units and their precision, the stage each strain is
  lost at). A hit is correct when its protein is at least ``--within`` identical to the
  strain's best indexed protein.
"""

import argparse
import contextlib
import gzip
import json
import math
import re
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tree-index"))
import anchor  # noqa: E402
import compact  # noqa: E402
import prepare  # noqa: E402

from kmer_functional_profiler import _core  # noqa: E402

TREE = prepare.ROOT / "workflows/tree-index"
DEPTHS = (0.5, 1, 2, 5, 10, 30)
IDENTITIES = ("1", "0.98", "0.95", "0.9", "0.85", "0.8", "0.75", "0.7")
SETS = ("strains", "relatives")
RELEASE_PROTEINS = 5.74e9
K, RATE = 9, 32  # defaults: layer k; one key per RATE representative residues, on average
RATES = (32, 8, 4)  # the index's sizes: standard, and the larger options (keys x4 and x8)
FLOOR = 30.0  # read-pair evidence (bits) before a candidate is pseudo-aligned
LOOSE = (0.7, 0.5)  # also scored: the representative within these identities of the best
MIN_READS = 0.5  # expected reads (EM weight) for a unit to be present
BTOP = re.compile(r"(\d+)|(..)")


def depth_runs(a: argparse.Namespace, name: str, targets: list[tuple[int, str]]) -> None:
    """One ``spike.py simulate`` per depth on its share of ``targets``, merged into
    ``<name>/reads_{1,2}.fq.gz`` and ``sim.parquet`` (strain key, depth, member, identity_bin);
    reads and strains renamed ``d<j>s<strain>``."""
    d = a.out / name
    d.mkdir(parents=True, exist_ok=True)
    rows = []
    with contextlib.ExitStack() as stack:
        outs = {m: stack.enter_context(gzip.open(d / f"reads_{m}.fq.gz", "wt")) for m in (1, 2)}
        for j, depth in enumerate(DEPTHS):
            simulate(a, name, d, j, depth, targets, outs, rows)
    pl.concat(rows).with_columns(member=pl.col("member").cast(pl.Int64)).write_parquet(
        d / "sim.parquet"
    )


def simulate(a: argparse.Namespace, name: str, d: Path, j: int, depth: float,
             targets: list, outs: dict, rows: list) -> None:  # fmt: skip
    """``spike.py simulate`` at ``depth`` exactly on every len(DEPTHS)-th target from j; its
    reads appended to ``outs`` and its truth to ``rows``, renamed ``d<j>``."""
    run = d / f"depth{j}"
    run.mkdir(exist_ok=True)
    share = targets[j :: len(DEPTHS)]
    (run / "targets.fa").write_text("".join(f">{i}\n{s}\n" for i, s in share))
    idents = IDENTITIES if name == "strains" else ("1",)
    prepare.run(prepare.PY, prepare.SPIKE, "simulate", "--targets", run / "targets.fa",
                "--present", 1, "--identities", *idents, "--substitution", "blosum62",
                "--depth-mu", math.log(depth), "--depth-sigma", 0, "--seed", a.seed + j,
                "--out", run)  # fmt: skip
    for m in (1, 2):
        with gzip.open(run / f"reads_{m}.fq.gz", "rt") as f:
            for i, line in enumerate(f):
                outs[m].write(f"@d{j}{line[1:]}" if i % 4 == 0 else line)
    t = pl.read_csv(run / "truth.tsv", separator="\t", columns=["strain", "member", "identity_bin"])
    rows.append(
        t.with_columns(strain=pl.format("d{}s{}", pl.lit(j), "strain"), depth=pl.lit(float(depth)))
    )
    shutil.rmtree(run)


def reads(a: argparse.Namespace) -> None:
    ref = pl.read_parquet(a.reference, columns=["protein_id", "cluster_rep", "sequence"])
    rng = np.random.default_rng(a.seed)
    pool = ref.filter(pl.col("sequence").str.len_chars() >= 60)
    pick = pool[rng.choice(pool.height, a.strains * len(DEPTHS), replace=False).tolist()]
    depth_runs(a, "strains", list(zip(pick["protein_id"], pick["sequence"], strict=True)))
    held = a.holdout / "targets.fa"
    recs = [r.split("\n", 1) for r in held.read_text().split(">")[1:]]
    recs = [(int(h), s.replace("\n", "")) for h, s in recs]
    depth_runs(a, "relatives", [recs[i] for i in rng.permutation(len(recs))])
    unit = ref.select(best="protein_id", unit="cluster_rep")
    rel = pl.read_csv(a.holdout / "relatives.tsv", separator="\t").select(
        member="protein_id", best="relative", aai="identity_bin"
    )
    for name in SETS:
        sim = pl.read_parquet(a.out / name / "sim.parquet")
        if name == "strains":
            truth = sim.with_columns(best="member", aai="identity_bin")
        else:
            truth = sim.join(rel, on="member")
        truth.select("strain", "depth", "aai", "member", "best").join(
            unit, on="best"
        ).write_parquet(a.out / name / "truth.parquet")


def tag(a: argparse.Namespace) -> str:
    """The index's name by its sparsity: ``r32`` (the default), ``r8`` or ``r4`` (the larger)."""
    return f"r{a.rate:g}"


def build(a: argparse.Namespace) -> None:
    """The compact index (``compact.py shard`` per shard, then ``merge``)."""
    w = a.out / f"index-{tag(a)}"
    times: dict[str, object] = {}
    py, cx = prepare.PY, TREE / "compact.py"
    for s in range(a.shards):
        if not (w / f"seqs/{s:04d}.json").exists():
            t0 = time.perf_counter()
            prepare.run(py, cx, "shard", "--members", a.reference, "--index", w, "--shards",
                        a.shards, "--shard", s, "--k", a.k, "--rate", a.rate)  # fmt: skip
            times[f"shard_{s}"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    prepare.run(py, cx, "merge", "--index", w, "--k", a.k, "--rate", a.rate)
    times["merge"] = time.perf_counter() - t0
    meta = json.loads((w / "meta.json").read_text())
    times |= {"postings": meta["postings"], "index_bytes": meta["bytes"],
              "seq_bytes": sum(f.stat().st_size for f in (w / "seqs").glob("*"))}  # fmt: skip
    (a.out / f"build_times_{tag(a)}.json").write_text(json.dumps(times, indent=1))
    print(json.dumps(times, indent=1))


def query(a: argparse.Namespace) -> None:
    scale = (
        a.chance_scale
        or RELEASE_PROTEINS / pl.scan_parquet(a.reference).select(pl.len()).collect().item()
    )
    times: dict[str, object] = {"chance_scale": scale}
    for name in SETS:
        d = a.out / name
        r = [d / "reads_1.fq.gz", d / "reads_2.fq.gz"]
        t0 = time.perf_counter()
        prepare.run(prepare.PY, TREE / "compact.py", "query", "--index", a.out / f"index-{tag(a)}",
                    "--reads", *r, "--out", d / f"c-{tag(a)}", "--floor", FLOOR,
                    "--chance-scale", scale)  # fmt: skip
        times[f"{name}_index_s"] = time.perf_counter() - t0  # the index's cost: tier 1, confirm
        t0 = time.perf_counter()
        prepare.run(prepare.PY, TREE / "recruit.py", "--reads", *r, "--targets",
                    d / f"c-{tag(a)}/confirmed.fa", "--seeded", d / f"c-{tag(a)}/seeded.tsv",
                    "--out", d / f"r-{tag(a)}",
                    "--chance-scale", scale)  # fmt: skip
        times[f"{name}_recruit_s"] = time.perf_counter() - t0  # a step of its own
        times[name] = json.loads((d / f"c-{tag(a)}/query.json").read_text())
        times[f"{name}_recruit"] = json.loads((d / f"r-{tag(a)}/recruit.json").read_text())
    (a.out / f"query_times_{tag(a)}.json").write_text(json.dumps(times, indent=1))
    print(json.dumps(times, indent=1))


def identities(pairs: pl.DataFrame, seqs: dict[int, str]) -> pl.DataFrame:
    """(a, b, identity) for each distinct protein pair: matches in the seed alignment of a on b
    (``_core.seed_align``) over the shorter protein's length (0 if none found)."""
    pairs = pairs.unique().sort("b", "a")
    ids_a, ids_b = pairs["a"].unique().sort().to_list(), pairs["b"].unique().sort().to_list()
    ia, ib = {p: i for i, p in enumerate(ids_a)}, {p: i for i, p in enumerate(ids_b)}
    pq = np.array([ia[x] for x in pairs["a"].to_list()], dtype=np.uint32)
    pt = np.array([ib[x] for x in pairs["b"].to_list()], dtype=np.uint32)
    found, *_, btop = _core.seed_align(
        [seqs[x].encode() for x in ids_a], [seqs[x].encode() for x in ids_b], pq, pt,
        anchor.blosum().astype(np.int8).ravel(), k=anchor.SEED_K, band=anchor.SEED_BAND,
        x_drop=anchor.X_DROP,
    )  # fmt: skip
    shorter = np.minimum(
        [len(seqs[x]) for x in pairs["a"].to_list()], [len(seqs[x]) for x in pairs["b"].to_list()]
    )
    match = np.array([sum(int(m[0]) for m in BTOP.findall(b) if m[0]) for b in btop])
    return pairs.with_columns(
        identity=pl.Series(np.where(found, match / np.maximum(shorter, 1), 0.0))
    )


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if not n:
        return (math.nan, math.nan)
    p, d = k / n, 1 + z * z / n
    c, h = (p + z * z / (2 * n)) / d, z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def score(a: argparse.Namespace) -> None:
    ref = pl.read_parquet(a.reference, columns=["protein_id", "cluster_rep", "sequence"])
    seqs = dict(ref.select("protein_id", "sequence").iter_rows())
    rep_of = dict(compact.unit_reps(ref).iter_rows())
    held = a.holdout / "targets.fa"
    seqs |= {
        int(h): s.replace("\n", "")
        for h, s in (r.split("\n", 1) for r in held.read_text().split(">")[1:])
    }
    recall_rows, target_rows = [], []
    for name in SETS:
        d = a.out / name
        truth = pl.read_parquet(d / "truth.parquet")
        names = anchor.read_names(d / "reads_1.fq.gz")
        read_strain = pl.DataFrame(
            {"read": np.arange(len(names)), "strain": [n.rsplit("_", 1)[0] for n in names]}
        )
        units = pl.read_csv(d / f"c-{tag(a)}/units.tsv", separator="\t")  # the confirmed
        # a unit is reported by its representative; each read pair's top units, ties kept
        top = (
            pl.read_parquet(d / f"r-{tag(a)}/reads.parquet", columns=["read", "target"])
            .select("read", protein_id=pl.col("target").cast(pl.Int64))
            .join(units.select("unit", "protein_id"), on="protein_id")
            .join(read_strain, on="read")
            .join(truth.select("strain", "best"), on="strain")
        )
        recruited = pl.read_csv(
            d / f"r-{tag(a)}/targets.tsv", separator="\t", schema_overrides={"target": pl.Int64}
        )
        present = (
            recruited.filter(pl.col("reads") >= MIN_READS)
            .select(protein_id="target")
            .join(units.select("unit", "protein_id"), on="protein_id")
        )
        # every present unit against the strains whose reads hit it
        support = present.join(top.select("unit", "best").unique(), on="unit", how="left")
        ident = identities(
            pl.concat([top.select(a="protein_id", b="best"),
                       support.drop_nulls("best").select(a="protein_id", b="best")]), seqs,
        )  # fmt: skip
        ok = ident.filter(pl.col("identity") >= a.within).select("a", "b", hit=pl.lit(True))
        top = top.join(ok, left_on=["protein_id", "best"], right_on=["a", "b"], how="left")
        support = support.join(ok, left_on=["protein_id", "best"], right_on=["a", "b"], how="left")
        per_read = top.group_by("read", "strain").agg(right=pl.col("hit").any().fill_null(False))
        found_best = support.filter(pl.col("hit").fill_null(False)).select("best").unique()
        # looser: a present unit's representative within each LOOSE identity of the best
        sup = support.drop_nulls("best").join(ident, left_on=["protein_id", "best"],
                                              right_on=["a", "b"], how="left")  # fmt: skip
        loose = {x: sup.filter(pl.col("identity") >= x)["best"].unique().implode() for x in LOOSE}
        conf, shown = set(units["protein_id"].to_list()), set(present["protein_id"].to_list())
        confirmed = [rep_of.get(b) in conf for b in truth["best"].to_list()]
        unit_found = [rep_of.get(b) in shown for b in truth["best"].to_list()]
        by_strain = per_read.group_by("strain").agg(assigned=pl.len(), right=pl.col("right").sum())
        strain_reads = read_strain.group_by("strain").agg(pairs=pl.len())
        st = (
            truth.with_columns(confirmed=pl.Series(confirmed, dtype=pl.Boolean),
                               unit_found=pl.Series(unit_found, dtype=pl.Boolean))
            .join(strain_reads, on="strain", how="left")
            .join(by_strain, on="strain", how="left")
            .with_columns(pl.col("pairs", "assigned", "right").fill_null(0),
                          found=pl.col("best").is_in(found_best["best"].implode()),
                          **{f"found_{x}": pl.col("best").is_in(v) for x, v in loose.items()})
        )  # fmt: skip
        st.write_parquet(d / f"route90_strains_{tag(a)}.parquet")
        for (depth, aai), g in st.group_by("depth", "aai"):
            n, k = g.height, int(g["found"].sum())
            lo, hi = wilson(k, n)
            recall_rows.append({
                "set": name, "depth": depth, "aai": aai, "strains": n, "target_recall": k / n,
                "lo": lo, "hi": hi, **{f"recall_{x}": float(g[f"found_{x}"].mean()) for x in LOOSE},
                "unit_recall": float(g["unit_found"].mean()),
                "confirmed": float(g["confirmed"].mean()),
                "pairs": int(g["pairs"].sum()),
                "read_recall": g["right"].sum() / max(g["pairs"].sum(), 1),
                "read_precision": g["right"].sum() / max(g["assigned"].sum(), 1),
            })  # fmt: skip
        true_pos = support.group_by("unit").agg(t=pl.col("hit").any().fill_null(False))
        target_rows.append({
            "set": name, "present": present.height, "true": int(true_pos["t"].sum()),
            "precision": float(true_pos["t"].mean()) if true_pos.height else math.nan,
            "unsupported": int(support.filter(pl.col("best").is_null()).height),
            "strains": truth.height, "found": int(st["found"].sum()),
            "lost_not_confirmed": int((~st["confirmed"]).sum()),
            "lost_after_confirm": int((st["confirmed"] & ~st["found"]).sum()),
        })  # fmt: skip
    rec = pl.DataFrame(recall_rows).sort("set", "depth", "aai", descending=[False, False, True])
    rec.write_csv(a.out / f"route90_recall_{tag(a)}.tsv", separator="\t", float_precision=3)
    pl.DataFrame(target_rows).write_csv(
        a.out / f"route90_targets_{tag(a)}.tsv", separator="\t", float_precision=3
    )
    with pl.Config(tbl_rows=200, tbl_cols=20):
        print(rec)
        print(pl.DataFrame(target_rows))


STEPS = {"reads": reads, "build": build, "query": query, "score": score}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("step", choices=[*STEPS, "all"])
    p.add_argument("--out", type=Path, default=prepare.ROOT / "data/route90")
    p.add_argument("--holdout", type=Path, default=prepare.ROOT / "data/pfam-anchor/holdout")
    p.add_argument("--reference", type=Path, help="default: <holdout>/members.parquet")
    p.add_argument("--strains", type=int, default=240, help="per depth")
    p.add_argument("--shards", type=int, default=4)
    p.add_argument(
        "--chance-scale", type=float, default=0.0, help="0: release / reference proteins"
    )
    p.add_argument(
        "--within", type=float, default=0.9, help="identity to the best target that counts"
    )
    p.add_argument("--k", type=int, default=K, help="tier 1 key length")
    p.add_argument(
        "--rate",
        type=float,
        default=RATE,
        choices=RATES,
        help="one key per RATE rep residues: 32, or 8 and 4 for the larger indexes",
    )
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--max-gb", type=float, default=prepare.MAX_GB)
    a = p.parse_args()
    a.reference = a.reference or a.holdout / "members.parquet"
    prepare.MAX_GB = a.max_gb
    prepare.watchdog(a.max_gb)
    for name, step in STEPS.items():
        if a.step in (name, "all"):
            step(a)


if __name__ == "__main__":
    main()
