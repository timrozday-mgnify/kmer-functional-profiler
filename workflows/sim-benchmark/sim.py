"""Simulated abundance benchmark: exact per-unit truth for the phase-4 models, run locally.

Reference: ``--families`` full-length seed proteins (200-600 aa) from ``--seed-fasta``,
by default the MGnify sample (``scripts/fetch_mgnify_sample.py``; ``FL=0`` records
skipped). Each family has ``--paralogs`` units: the seed mutated to 85-95% identity, so
paralogs share k-mers as close MGnify90 clusters do. A unit's
``--members`` members are its centroid mutated to 97%, as in a 90% cluster.

Sample (one per seed): ``--present`` of the units, each carried by one strain, the unit
centroid mutated to an identity drawn from ``IDENTITIES``; depth is lognormal. Strains are
back-translated with random codons between random flanks; single-end reads of ``READ`` bp
from either strand get ``--error`` substitutions, plus ``--decoys`` random-DNA reads.

Every index config (``CONFIGS``) is queried with every sample, and every assignment rule
(count and abundance pair, as in the fmh benchmark) is scored at ``min_hits`` 1 against
the true depth:

- ``purity``, ``completeness``: of units detected;
- ``spearman_tp``: rank correlation of estimate and depth over true positives;
- ``l1``, ``l1_family``: L1 between relative abundances over units, and over families
  (sums per family: what a Pfam-level report sees);
- ``log_ratio_sd``: sd of log(estimate / depth) over true positives; 0 = one scale factor;
- ``bias_<identity>``: median estimate / depth of true positives at that strain identity,
  over the median at 100%; 1 = divergence costs nothing.

Phase 7 adds, per ``--frames`` mode (one set of rows each, column ``frames``):

- ``complete_<bin>``, ``bias_len_<bin>``: completeness, and median estimate / depth of true
  positives over the median of all of them, by protein length (``LENGTH_BINS`` aa). Reads
  crossing a gene end are dropped by the stop-free filter, which costs short proteins more;
  ``--frames edges`` is meant to recover them.
- On ``gather_zi`` rows, containment AAI against the truth (substitutions only, so identity
  is per position): ``aai_bias_<identity>`` (median ``aai`` minus the strain's mean identity
  to its unit's members, by the strain's identity to the unit centroid);
  ``aai_bias_nearest`` and ``aai_bias_centroid`` (against the nearest member and the
  centroid, all identities); ``aai_spearman``; ``aai_cover`` and ``aai_cover_centroid``
  (the posterior interval holds the mean-member and the centroid identity; with
  ``--draws``); ``aai_naive_bias``; ``aai_n`` (true positives with an ``aai``).

``index_mb`` is the size of the index's lookup tables. The ``gather_zi`` rows add how
``present_prob`` separates true from false positives (:func:`presence_scores`) and, with
``--draws`` D, interval calibration (:func:`calibration`). Writes ``scores.tsv`` (per
seed) and ``summary.tsv`` (means) to ``--out`` and prints the summary. Seconds per seed.
"""

import argparse
import itertools
import math
import random
from pathlib import Path

import numpy as np
import polars as pl

from kmer_functional_profiler.index import Index, IndexParams, build_index
from kmer_functional_profiler.query import profile
from kmer_functional_profiler.reference import BASES, CODE_11, reverse_complement

ROOT = Path(__file__).resolve().parents[2]
SEEDS_FASTA = ROOT / "data" / "mgnify" / "mgy_clusters_head.faa"
AMINO = "ACDEFGHIKLMNPQRSTVWY"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
SYNONYMS = {aa: [c for c, t in zip(CODONS, CODE_11, strict=True) if t == aa] for aa in AMINO}
IDENTITIES = (1.0, 0.95, 0.9, 0.85)
LENGTH_BINS = ((200, 300), (300, 450), (450, 601))  # aa, half-open
READ, FLANK = 150, 150
CONFIGS = {
    "dense": IndexParams(t_base=1.0, n_min=0),
    "s10": IndexParams(t_base=0.1, n_min=0),
    "floor": IndexParams(),  # MGnify defaults: t_base 0.001, n_min 8, t_cap 0.2
    # plus a dense tier at 1 in 50, 20, 5 and 1
    **{f"floor_d{round(t * 100)}": IndexParams(t_dense=t) for t in (0.02, 0.05, 0.2, 1.0)},
}
# Detection count and its abundance estimate, per assignment rule.
RULES = {
    "none": ("kmers_hit", "coverage"),
    "gather_em": ("kmers_unique", "coverage_em"),
    "gather_zi": ("kmers_unique", "coverage_zi"),
    "gather_zib": ("kmers_unique", "coverage_zib"),
    "gather_zip": ("kmers_unique", "coverage_zip"),
    "gather_zi_copies": ("kmers_unique", "abundance_zi"),
    "wta": ("kmers_wta", "coverage_wta"),
    "ufirst": ("kmers_ufirst", "coverage_ufirst"),
}


def seed_proteins(path: Path, n: int) -> list[str]:
    records = path.read_text().split(">")[1:]
    full = []
    for r in records:
        header, seq = r.split("\n", 1)
        seq = seq.replace("\n", "")
        if "FL=0" not in header and 200 <= len(seq) <= 600 and set(seq) <= set(AMINO):
            full.append(seq)
    return full[:n]


def mutate(seq: str, identity: float, rng: random.Random) -> str:
    out = list(seq)
    for i in rng.sample(range(len(seq)), round((1 - identity) * len(seq))):
        out[i] = rng.choice(AMINO.replace(seq[i], ""))
    return "".join(out)


def reference(args: argparse.Namespace) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Members table for the index, and the unit table (``unit``, ``family``, ``centroid``,
    ``twin_of``: the unit a near-identical twin copies, else null)."""
    rng = random.Random(0)
    units: list[tuple[int, int, str, int | None]] = []
    members: list[tuple[int, int, bool, str]] = []

    def add(family: int, centroid: str, twin_of: int | None) -> None:
        unit = len(units)
        units.append((unit, family, centroid, twin_of))
        first_id = len(members)
        members.extend(
            (first_id + i, unit, True, mutate(centroid, 0.97, rng)) for i in range(args.members)
        )

    for family, seed in enumerate(seed_proteins(args.seed_fasta, args.families)):
        first = len(units)
        for _ in range(args.paralogs):
            add(family, mutate(seed, rng.uniform(0.85, 0.95), rng), None)
        if rng.random() < args.twins:  # a near-identical paralog: barely separable
            add(family, mutate(units[first][2], 0.99, rng), first)
    return (
        pl.DataFrame(
            members, schema=["protein_id", "cluster_rep", "full_length", "sequence"], orient="row"
        ),
        pl.DataFrame(
            units,
            schema={
                "unit": pl.Int64,
                "family": pl.Int64,
                "centroid": pl.String,
                "twin_of": pl.Int64,
            },
            orient="row",
        ),
    )


def sample(units: pl.DataFrame, args: argparse.Namespace, seed: int, fasta: Path) -> pl.DataFrame:
    """Write the sample's reads to ``fasta``; return truth (``unit``, ``depth``, ``identity``,
    ``strain``)."""
    rng = random.Random(seed)
    chosen = rng.sample(range(units.height), round(args.present * units.height))
    truth, reads = [], []
    for unit in chosen:
        identity = rng.choice(IDENTITIES)
        depth = min(math.exp(rng.gauss(1.0, 1.0)), 50.0)
        strain = mutate(units["centroid"][unit], identity, rng)
        cds = "".join(rng.choice(SYNONYMS[aa]) for aa in strain) + "TAA"
        seq = "".join(rng.choices(BASES, k=FLANK)) + cds + "".join(rng.choices(BASES, k=FLANK))
        n = round(depth * len(seq) / READ)
        for _ in range(n):
            start = rng.randrange(len(seq) - READ + 1)
            read = seq[start : start + READ]
            if rng.random() < 0.5:
                read = reverse_complement(read.encode()).decode()
            reads.append(read)
        truth.append((unit, depth, identity, strain))
    n_decoys = round(args.decoys * len(reads))
    reads += ["".join(rng.choices(BASES, k=READ)) for _ in range(n_decoys)]
    reads = [
        "".join(rng.choice(BASES.replace(b, "")) if rng.random() < args.error else b for b in r)
        for r in reads
    ]
    fasta.write_text("".join(f">r{i}\n{r}\n" for i, r in enumerate(reads)))
    return pl.DataFrame(truth, schema=["unit", "depth", "identity", "strain"], orient="row")


def score(truth: pl.DataFrame, found: pl.DataFrame, families: pl.DataFrame) -> dict[str, float]:
    """Scores of one rule; ``found`` has ``unit`` and ``estimate`` for the detected units;
    ``families`` has ``unit``, ``family`` and ``length`` (aa)."""
    both = (
        truth.select("unit", "depth", "identity")
        .join(found, on="unit", how="full", coalesce=True)
        .join(families, on="unit")
        .with_columns(pl.col("depth", "estimate").fill_null(0.0))
    )
    tp = both.filter(pl.col("depth") > 0, pl.col("estimate") > 0).with_columns(
        ratio=pl.col("estimate") / pl.col("depth")
    )

    def l1(df: pl.DataFrame) -> float:
        return float(
            (df["depth"] / df["depth"].sum() - df["estimate"] / max(df["estimate"].sum(), 1e-300))
            .abs()
            .sum()
        )

    ref = tp.filter(pl.col("identity") == 1.0)["ratio"].median()
    return {
        "purity": tp.height / max(found.height, 1),
        "completeness": tp.height / truth.height,
        "spearman_tp": tp.select(pl.corr("depth", "estimate", method="spearman")).item(),
        "l1": l1(both),
        "l1_family": l1(both.group_by("family").agg(pl.col("depth", "estimate").sum())),
        "log_ratio_sd": float(np.log(tp["ratio"].to_numpy()).std()),
        **{
            f"bias_{i}": tp.filter(pl.col("identity") == i)["ratio"].median() / ref  # type: ignore[operator]
            for i in IDENTITIES[1:]
        },
        **length_scores(both, tp),
    }


def length_scores(both: pl.DataFrame, tp: pl.DataFrame) -> dict[str, float | None]:
    """Completeness and relative abundance bias of true positives by protein length."""
    out: dict[str, float | None] = {}
    overall = tp["ratio"].median()
    for lo, hi in LENGTH_BINS:
        in_bin = pl.col("length").is_between(lo, hi, closed="left")
        present = both.filter(pl.col("depth") > 0, in_bin)
        found = tp.filter(in_bin)
        out[f"complete_{lo}"] = found.height / present.height if present.height else None
        out[f"bias_len_{lo}"] = (
            found["ratio"].median() / overall if found.height and overall else None  # type: ignore[operator]
        )
    return out


def _identity(a: str, b: str) -> float:
    """Share of equal positions (sequences of one length: substitutions only)."""
    return sum(x == y for x, y in zip(a, b, strict=True)) / len(a)


def aai_scores(
    truth: pl.DataFrame, result: pl.DataFrame, members: pl.DataFrame
) -> dict[str, float | None]:
    """Containment AAI of the true positives against their strains' identity to the unit's
    members (mean and nearest), by the strain's identity to the unit centroid."""
    if "aai" not in result.columns:
        return {}
    by_unit = members.group_by(unit=pl.col("cluster_rep").cast(pl.Int64)).agg("sequence")
    rows = []
    joined = truth.join(result.filter(pl.col("aai").is_not_null()), on="unit").join(
        by_unit, on="unit"
    )
    for row in joined.iter_rows(named=True):
        ids = [_identity(row["strain"], m) for m in row["sequence"]]
        rows.append(
            {
                "identity": row["identity"],
                "aai": row["aai"],
                "naive": row["aai_naive"],
                "mean": sum(ids) / len(ids),
                "nearest": max(ids),
                "lo": row.get("aai_lo"),
                "hi": row.get("aai_hi"),
            }
        )
    if not rows:
        return {"aai_n": 0}
    df = pl.DataFrame(rows, infer_schema_length=None)
    out: dict[str, float | None] = {
        "aai_n": df.height,
        "aai_bias_nearest": (df["aai"] - df["nearest"]).median(),  # type: ignore[dict-item]
        "aai_bias_centroid": (df["aai"] - df["identity"]).median(),  # type: ignore[dict-item]
        "aai_naive_bias": (df["naive"] - df["mean"]).median(),  # type: ignore[dict-item]
        "aai_spearman": df.select(pl.corr("aai", "mean", method="spearman")).item()
        if df.height > 2
        else None,
        **{
            f"aai_cover{suffix}": df.drop_nulls("lo")
            .select(pl.col(truth_col).is_between("lo", "hi").mean())
            .item()
            if df["lo"].drop_nulls().len()
            else None
            for suffix, truth_col in (("", "mean"), ("_centroid", "identity"))
        },
    }
    for i in IDENTITIES:
        part = df.filter(pl.col("identity") == i)
        out[f"aai_bias_{i}"] = (part["aai"] - part["mean"]).median() if part.height else None  # type: ignore[assignment]
    return out


def presence_scores(truth: pl.DataFrame, result: pl.DataFrame) -> dict[str, float]:
    """``present_prob`` of the units gather keeps: mean over true (``prob_tp``) and false
    (``prob_fp``) positives, and the share of each below 0.5 (``flag_tp``, ``flag_fp``)."""
    d = result.filter(pl.col("kmers_unique") >= 1).select(
        "present_prob", tp=pl.col("unit").is_in(truth["unit"].implode())
    )
    out = {}
    for name, part in (("tp", d.filter("tp")), ("fp", d.filter(~pl.col("tp")))):
        out[f"prob_{name}"] = part["present_prob"].mean()
        out[f"flag_{name}"] = (part["present_prob"] < 0.5).mean()
    return out  # type: ignore[return-value]


def calibration(truth: pl.DataFrame, result: pl.DataFrame, units: pl.DataFrame) -> dict[str, float]:
    """How well ``coverage_zi``'s posterior intervals and ambiguity groups hold the truth.

    Truth is read depth; k-mer coverage is depth times a constant (read length, k, errors),
    taken as the median estimate / depth over well-covered strains at 100% identity.
    ``ci_cover`` is over true positives, ``ci_cover_low`` over those at depth <= 2,
    ``ci_cover_twin`` over those with a detected near-identical twin, and ``ci_width`` the
    median log(hi / lo). ``group_cover`` is how often a group's interval holds its members'
    true total; ``grouped`` the share of true positives in a group; ``twins_grouped`` the
    share of detected twin pairs placed in one group.
    """
    detected = result.filter(pl.col("kmers_unique") >= 1)
    tp = truth.join(detected, on="unit").filter(pl.col("coverage_zi") > 0)
    scale = (
        tp.filter(pl.col("identity") == 1.0, pl.col("depth") > 5)
        .select((pl.col("coverage_zi") / pl.col("depth")).median())
        .item()
    )
    pairs = units.filter(pl.col("twin_of").is_not_null()).select(a="twin_of", b="unit")
    found = detected["unit"].to_list()
    pairs = pairs.filter(pl.col("a").is_in(found), pl.col("b").is_in(found))
    twinned = pairs["a"].to_list() + pairs["b"].to_list()
    group = dict(detected.select("unit", "ambiguity_group").iter_rows())
    groups = (
        detected.filter(pl.col("ambiguity_group").is_not_null())
        .join(truth, on="unit", how="left")
        .group_by("ambiguity_group")
        .agg(
            total=pl.col("depth").fill_null(0).sum() * scale,
            lo=pl.col("group_coverage_zi_lo").first(),
            hi=pl.col("group_coverage_zi_hi").first(),
        )
    )
    held = tp.select(
        twin=pl.col("unit").is_in(twinned),
        grouped=pl.col("ambiguity_group").is_not_null(),
        low=pl.col("depth") <= 2,
        inside=(pl.col("coverage_zi_lo") <= pl.col("depth") * scale)
        & (pl.col("depth") * scale <= pl.col("coverage_zi_hi")),
        width=(pl.col("coverage_zi_hi") / pl.col("coverage_zi_lo")).log(),
    )
    return {
        "ci_cover": held["inside"].mean(),  # type: ignore[dict-item]
        "ci_cover_low": held.filter("low")["inside"].mean(),  # type: ignore[dict-item]
        "ci_cover_twin": held.filter("twin")["inside"].mean(),  # type: ignore[dict-item]
        "ci_width": held["width"].median(),  # type: ignore[dict-item]
        "grouped": held["grouped"].mean(),  # type: ignore[dict-item]
        "group_cover": groups.select(pl.col("total").is_between("lo", "hi").mean()).item(),
        "twins_grouped": sum(
            group[a] is not None and group[a] == group[b] for a, b in pairs.iter_rows()
        )
        / max(pairs.height, 1),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-fasta", type=Path, default=SEEDS_FASTA)
    parser.add_argument("--families", type=int, default=100)
    parser.add_argument("--paralogs", type=int, default=3)
    parser.add_argument("--members", type=int, default=4)
    parser.add_argument("--present", type=float, default=0.5)
    parser.add_argument("--error", type=float, default=0.002)
    parser.add_argument("--decoys", type=float, default=0.2)
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--draws", type=int, default=0, help="posterior draws for intervals")
    parser.add_argument("--twins", type=float, default=0.0, help="families with a 99%% twin")
    parser.add_argument("--configs", nargs="+", default=list(CONFIGS), choices=list(CONFIGS))
    parser.add_argument("--frames", nargs="+", default=["stopfree"],
                        help="frame modes queried, e.g. stopfree edges:20 all")  # fmt: skip
    parser.add_argument("--out", type=Path, default=Path("sim-results"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    members, units = reference(args)
    members.write_parquet(args.out / "members.parquet")
    families = units.select("unit", "family", length=pl.col("centroid").str.len_chars())
    rows = []
    for config in args.configs:
        build_index(args.out / "members.parquet", args.out / f"index_{config}", CONFIGS[config])
        index = Index.load(args.out / f"index_{config}")
        stats = index.meta["stats"]
        index_mb = (stats["tier2_bytes"] + stats.get("dense_bytes", 0)) / 1e6
        for seed, frames in itertools.product(range(1, args.seeds + 1), args.frames):
            reads = args.out / f"reads_{seed}.fa"
            truth = sample(units, args, seed, reads)
            result = profile(
                index, reads, frames=frames, draws=args.draws, all_estimators=True
            ).with_columns(unit=pl.col("cluster_rep").cast(pl.Int64))
            for rule, (count, abundance) in RULES.items():
                found = result.filter(pl.col(count) >= 1).select("unit", estimate=abundance)
                extra = (
                    {
                        **presence_scores(truth, result),
                        **(calibration(truth, result, units) if args.draws else {}),
                        **aai_scores(truth, result, members),
                    }
                    if rule == "gather_zi"
                    else {}
                )
                key = {"config": config, "frames": frames, "rule": rule, "seed": seed}
                rows.append({**key, "index_mb": index_mb, **score(truth, found, families), **extra})
    scores = pl.DataFrame(rows, infer_schema_length=None)
    scores.write_csv(args.out / "scores.tsv", separator="\t")
    summary = (
        scores.group_by("config", "frames", "rule", maintain_order=True)
        .agg(pl.exclude("seed").mean())
        .with_columns(pl.selectors.float().round(3))
    )
    summary.write_csv(args.out / "summary.tsv", separator="\t")
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=200):
        print(summary)


if __name__ == "__main__":
    main()
