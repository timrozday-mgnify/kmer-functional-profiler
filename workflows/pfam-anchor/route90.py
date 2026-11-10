"""The tree index on its own, with no anchoring: tier 1 on MGnify90 units (~90% AAI), each
unit's keys on one representative (its medoid), k 9, about 1 in 8 of the representatives'
residues keyed and spent per unit in proportion to sqrt(members), a read-pair floor of 30
bits; then each candidate (read pair, unit) pseudo-aligned to the representative to confirm
it (``verify.py``), every leaf of the confirmed units pseudo-aligned (``tier2.py``'s dense
routing, then its anchored diagonal-band scan only: no gapped stage, no neighbouring or
merely routed units), and the top hits fully aligned (``_core.seed_align``, ``call``'s
DIAMOND lines). Recall and precision by depth x AAI (docs/tree-index-plan.md, Progress log:
the route on its own, 2026-11-09).

Steps, into ``--out``:

- ``reads``: ``strains/`` (reference proteins mutated to 100-70% identity, BLOSUM62) and
  ``relatives/`` (pfam-anchor's held-out proteins, at their real identity to the nearest
  protein left in the reference), each simulated at every depth of :data:`DEPTHS` (one
  ``spike.py simulate`` per depth, at that depth exactly) and merged; ``truth.parquet`` per set:
  strain, depth, aai, best (the indexed protein a perfect method would report).
- ``build``: the index on pfam-anchor's reference (``partition.py``, ``build.py`` per part,
  ``merge.py``, ``store.py``).
- ``query``: per set, ``query.py --decoy``, ``verify.py --floor 30``, ``tier2.py``, then the
  seed alignments of each read pair's top leaves (``seed_{1,2}.tsv``); bars at ``--chance-scale``
  (the release's proteins over the reference's, by default).
- ``score``: ``route90_recall.tsv`` (per set x depth x AAI: target recall with Wilson intervals
  over strains, read recall and precision), ``route90_targets.tsv`` (present targets and their
  precision, the stage each strain is lost at) and ``route90_cost.json``. A hit is correct when
  its protein is at least ``--within`` identical to the strain's best indexed protein.
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
import anchor  # noqa: E402
import prepare  # noqa: E402

from kmer_functional_profiler import _core  # noqa: E402

TREE = prepare.ROOT / "workflows/tree-index"
DEPTHS = (0.5, 1, 2, 5, 10, 30)
IDENTITIES = ("1", "0.98", "0.95", "0.9", "0.85", "0.8", "0.75", "0.7")
SETS = ("strains", "relatives")
RELEASE_PROTEINS = 5.74e9
K, RATE = 9, 8  # defaults: layer k; one key per RATE representative residues, on average
POOL_WINDOW = 4  # the layer's window minimizers (~2/(w+1) of k-mers), thinned by sqrt:<mean>
FLOOR = 30.0  # read-pair evidence (bits) before a candidate is pseudo-aligned
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


def build(a: argparse.Namespace) -> None:
    w = a.out / "index"
    ref = pl.read_parquet(a.reference, columns=["protein_id", "cluster_rep", "sequence"])
    reps = ref.filter(pl.col("protein_id") == pl.col("cluster_rep"))
    mean = float((reps["sequence"].str.len_chars() - a.k + 1).clip(1).mean()) / a.rate
    times = {"keys_per_unit": mean}

    def step(key: str, *cmd: object) -> None:
        t0 = time.perf_counter()
        prepare.run(*cmd)
        times[key] = time.perf_counter() - t0

    py = prepare.PY
    if not (w / "parts/partition.json").exists():
        step("partition", py, TREE / "partition.py", "--members", a.reference, "--out", w / "parts",
             "--part-proteins", a.part_proteins)  # fmt: skip
    parts = sorted((w / "parts").glob("part*.parquet"))
    layer = ["--mode", "standard", "--layers", f"protein:{a.k}:{POOL_WINDOW}:rare",
             "--layer-groups", "units", "--layer-rep", "medoid", "--layer-keys", f"sqrt:{mean:.3f}",
             "--t-base", 0, "--window", 0, "--ladder=", "--n-1", 0]  # fmt: skip
    for part in parts:
        if not (w / "build" / part.stem / "proteins.parquet").exists():
            step(f"build_{part.stem}", py, TREE / "build.py", "--members", part, "--out",
                 w / "build" / part.stem, "--database-proteins", ref.height, *layer)  # fmt: skip
    if not (w / "index/parts.json").exists():
        step("merge", py, TREE / "merge.py", "--parts", *sorted((w / "build").glob("part*")),
             "--out", w / "index")  # fmt: skip
    if not (w / "store").exists():
        for part in parts:
            step(f"store_{part.stem}", py, TREE / "store.py", "build",
                 "--proteins", w / "build" / part.stem / "proteins.parquet", "--members", part,
                 "--parts", w / "index/parts.json", "--name", part.stem, "--format", "flat",
                 "--out", w / "store-parts" / part.stem)  # fmt: skip
        (w / "store").mkdir()
        for f in (w / "store-parts").glob("*/*"):
            shutil.copy(f, w / "store" / f.name)
    sizes = {f.name: f.stat().st_size for f in (w / "index").glob("*.npy")}
    times |= {
        "index_bytes": sizes,
        "postings": int(np.load(w / "index/hash.npy", mmap_mode="r").shape[0]),
    }
    (a.out / "build_times.json").write_text(json.dumps(times, indent=1))
    print(json.dumps({k: v for k, v in times.items() if k != "index_bytes"}, indent=1))


def best_leaves(index: Path, t2: Path) -> pl.DataFrame:
    """Each read pair's top leaves (``tier2.py``'s ``reads.parquet``: ties within its ambiguity)
    as (read, pos, protein_id), a leaf being one protein (:func:`leaf_protein`)."""
    return pl.read_parquet(t2 / "reads.parquet", columns=["read", "pos"]).join(
        leaf_protein(index), on="pos"
    )


def leaf_protein(index: Path) -> pl.DataFrame:
    # ponytail: a position's first protein (equal sampled sets: near-identical); the longest if
    # fragments start to matter
    prot = pl.read_parquet(index / "proteins.parquet", columns=["pos", "protein_id"])
    return prot.filter(pl.col("pos") >= 0).sort("pos", "protein_id").unique("pos", keep="first")


def query(a: argparse.Namespace) -> None:
    w = a.out / "index"
    idx, store = w / "index", w / "store"
    scale = (
        a.chance_scale
        or RELEASE_PROTEINS / pl.scan_parquet(a.reference).select(pl.len()).collect().item()
    )
    times = {"chance_scale": scale}
    py = prepare.PY
    seqs = dict(pl.read_parquet(a.reference, columns=["protein_id", "sequence"]).iter_rows())
    for name in SETS:
        d = a.out / name
        r = [d / "reads_1.fq.gz", d / "reads_2.fq.gz"]
        q, t2 = d / "q", d / "t2"
        steps = [
            ("tier1", [py, TREE / "query.py", "--index", idx, "--reads", *r, "--decoy",
                       "--chance-scale", scale, "--out", q]),
            ("confirm", [py, TREE / "verify.py", "--index", idx, "--store", store, "--reads", *r,
                         "--query", q, "--floor", FLOOR, "--chance-scale", scale]),
            ("leaves", [py, TREE / "tier2.py", "--index", idx, "--store", store, "--reads", *r,
                        "--query", q, "--chance-scale", scale, "--memory-gb", 3, "--no-handoff",
                        "--method", "anchored", "--neighbours", 0, "--route-below", 0,
                        "--out", t2]),
        ]  # fmt: skip
        for key, cmd in steps:
            t0 = time.perf_counter()
            prepare.run(*cmd)
            times[f"{name}_{key}_s"] = time.perf_counter() - t0
        names = anchor.read_names(r[0])
        mates = {}
        for m in (1, 2):
            with gzip.open(r[m - 1], "rt") as f:
                mates[m] = [line.strip() for i, line in enumerate(f) if i % 4 == 1]
        routed = (
            best_leaves(idx, t2).select("read", "protein_id").unique().sort("read", "protein_id")
        )
        t0 = time.perf_counter()
        out = anchor.rust_seed_lines(
            routed, mates, names, seqs, {p: p for p in seqs}, len(mates[1][0])
        )
        times[f"{name}_align_s"] = time.perf_counter() - t0
        times[f"{name}_lines"] = len(out[1]) + len(out[2])
        times[f"{name}_pairs"] = len(names)
        for m in (1, 2):
            (t2 / f"seed_{m}.tsv").write_text("".join(out[m]))
    (a.out / "query_times.json").write_text(json.dumps(times, indent=1))
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
    idx = a.out / "index/index"
    seqs = dict(pl.read_parquet(a.reference, columns=["protein_id", "sequence"]).iter_rows())
    held = a.holdout / "targets.fa"
    seqs |= {
        int(h): s.replace("\n", "")
        for h, s in (r.split("\n", 1) for r in held.read_text().split(">")[1:])
    }
    leaf = leaf_protein(idx)
    group = np.load(idx / "layer_group.npy").astype(np.int64)
    pos_of = dict(
        pl.read_parquet(idx / "proteins.parquet")
        .filter(pl.col("pos") >= 0)
        .select("protein_id", "pos")
        .iter_rows()
    )
    recall_rows, target_rows = [], []
    for name in SETS:
        d = a.out / name
        truth = pl.read_parquet(d / "truth.parquet")
        names = anchor.read_names(d / "reads_1.fq.gz")
        read_strain = pl.DataFrame(
            {"read": np.arange(len(names)), "strain": [n.rsplit("_", 1)[0] for n in names]}
        )
        top = (
            best_leaves(idx, d / "t2")
            .join(read_strain, on="read")
            .join(truth.select("strain", "best"), on="strain")
        )
        present = pl.DataFrame(
            {"pos": np.load(d / "t2/cand_verify_0.01.npy").astype(np.int64)}
        ).join(leaf, on="pos")
        # every present leaf against the strains whose reads hit it
        support = present.join(top.select("pos", "best").unique(), on="pos", how="left")
        ident = identities(
            pl.concat([top.select(a="protein_id", b="best"),
                       support.drop_nulls("best").select(a="protein_id", b="best")]), seqs,
        )  # fmt: skip
        ok = ident.filter(pl.col("identity") >= a.within).select("a", "b", hit=pl.lit(True))
        top = top.join(ok, left_on=["protein_id", "best"], right_on=["a", "b"], how="left")
        support = support.join(ok, left_on=["protein_id", "best"], right_on=["a", "b"], how="left")
        per_read = top.group_by("read", "strain").agg(right=pl.col("hit").any().fill_null(False))
        found_best = support.filter(pl.col("hit").fill_null(False)).select("best").unique()
        conf = set(np.load(d / "q/cand_verify_0.01.npy").astype(np.int64).tolist())
        confirmed = [pos_of.get(b, -1) >= 0 and int(group[pos_of[b]]) in conf
                     for b in truth["best"].to_list()]  # fmt: skip
        by_strain = per_read.group_by("strain").agg(assigned=pl.len(), right=pl.col("right").sum())
        strain_reads = read_strain.group_by("strain").agg(pairs=pl.len())
        st = (
            truth.with_columns(confirmed=pl.Series(confirmed, dtype=pl.Boolean))
            .join(strain_reads, on="strain", how="left")
            .join(by_strain, on="strain", how="left")
            .with_columns(pl.col("pairs", "assigned", "right").fill_null(0),
                          found=pl.col("best").is_in(found_best["best"].implode()))
        )  # fmt: skip
        st.write_parquet(d / "route90_strains.parquet")
        for (depth, aai), g in st.group_by("depth", "aai"):
            n, k = g.height, int(g["found"].sum())
            lo, hi = wilson(k, n)
            recall_rows.append({
                "set": name, "depth": depth, "aai": aai, "strains": n, "target_recall": k / n,
                "lo": lo, "hi": hi, "confirmed": float(g["confirmed"].mean()),
                "pairs": int(g["pairs"].sum()),
                "read_recall": g["right"].sum() / max(g["pairs"].sum(), 1),
                "read_precision": g["right"].sum() / max(g["assigned"].sum(), 1),
            })  # fmt: skip
        true_pos = support.group_by("pos").agg(t=pl.col("hit").any().fill_null(False))
        target_rows.append({
            "set": name, "present": present.height, "true": int(true_pos["t"].sum()),
            "precision": float(true_pos["t"].mean()) if true_pos.height else math.nan,
            "unsupported": int(support.filter(pl.col("best").is_null()).height),
            "strains": truth.height, "found": int(st["found"].sum()),
            "lost_not_confirmed": int((~st["confirmed"]).sum()),
            "lost_after_confirm": int((st["confirmed"] & ~st["found"]).sum()),
        })  # fmt: skip
    rec = pl.DataFrame(recall_rows).sort("set", "depth", "aai", descending=[False, False, True])
    rec.write_csv(a.out / "route90_recall.tsv", separator="\t", float_precision=3)
    pl.DataFrame(target_rows).write_csv(
        a.out / "route90_targets.tsv", separator="\t", float_precision=3
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
    p.add_argument("--part-proteins", type=int, default=40_000)
    p.add_argument(
        "--chance-scale", type=float, default=0.0, help="0: release / reference proteins"
    )
    p.add_argument(
        "--within", type=float, default=0.9, help="identity to the best target that counts"
    )
    p.add_argument("--k", type=int, default=K, help="tier 1 key length")
    p.add_argument("--rate", type=float, default=RATE, help="one key per RATE rep residues")
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
