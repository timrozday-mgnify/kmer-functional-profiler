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

``index_mb`` is the size of the index's lookup tables. With ``--bootstrap`` B, the
``gather_zi`` rows add interval calibration (:func:`calibration`). Writes ``scores.tsv`` (per
seed) and ``summary.tsv`` (means) to ``--out`` and prints the summary. Seconds per seed.
"""

import argparse
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
    """Members table for the index, and the unit table (``unit``, ``family``, ``centroid``)."""
    rng = random.Random(0)
    units, members = [], []
    for family, seed in enumerate(seed_proteins(args.seed_fasta, args.families)):
        for _ in range(args.paralogs):
            centroid = mutate(seed, rng.uniform(0.85, 0.95), rng)
            unit = len(units)
            units.append((unit, family, centroid))
            members += [
                (len(members), unit, True, mutate(centroid, 0.97, rng)) for _ in range(args.members)
            ]
    return (
        pl.DataFrame(
            members, schema=["protein_id", "cluster_rep", "full_length", "sequence"], orient="row"
        ),
        pl.DataFrame(units, schema=["unit", "family", "centroid"], orient="row"),
    )


def sample(units: pl.DataFrame, args: argparse.Namespace, seed: int, fasta: Path) -> pl.DataFrame:
    """Write the sample's reads to ``fasta``; return truth (``unit``, ``depth``, ``identity``)."""
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
        truth.append((unit, depth, identity))
    n_decoys = round(args.decoys * len(reads))
    reads += ["".join(rng.choices(BASES, k=READ)) for _ in range(n_decoys)]
    reads = [
        "".join(rng.choice(BASES.replace(b, "")) if rng.random() < args.error else b for b in r)
        for r in reads
    ]
    fasta.write_text("".join(f">r{i}\n{r}\n" for i, r in enumerate(reads)))
    return pl.DataFrame(truth, schema=["unit", "depth", "identity"], orient="row")


def score(truth: pl.DataFrame, found: pl.DataFrame, families: pl.DataFrame) -> dict[str, float]:
    """Scores of one rule; ``found`` has ``unit`` and ``estimate`` for the detected units."""
    both = (
        truth.join(found, on="unit", how="full", coalesce=True)
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
    }


def calibration(truth: pl.DataFrame, result: pl.DataFrame) -> dict[str, float]:
    """How often ``coverage_zi``'s bootstrap interval holds the true k-mer coverage.

    Truth is read depth; k-mer coverage is depth times a constant (read length, k, errors),
    taken as the median estimate / depth over well-covered strains at 100% identity.
    ``ci_cover`` is over true positives, ``ci_cover_low`` over those at depth <= 2, and
    ``ci_width`` the median log(hi / lo).
    """
    tp = truth.join(result.filter(pl.col("kmers_unique") >= 1), on="unit").filter(
        pl.col("coverage_zi") > 0
    )
    scale = (
        tp.filter(pl.col("identity") == 1.0, pl.col("depth") > 5)
        .select((pl.col("coverage_zi") / pl.col("depth")).median())
        .item()
    )
    held = tp.select(
        low=pl.col("depth") <= 2,
        inside=(pl.col("coverage_zi_lo") <= pl.col("depth") * scale)
        & (pl.col("depth") * scale <= pl.col("coverage_zi_hi")),
        width=(pl.col("coverage_zi_hi") / pl.col("coverage_zi_lo")).log(),
    )
    return {
        "ci_cover": held["inside"].mean(),  # type: ignore[dict-item]
        "ci_cover_low": held.filter("low")["inside"].mean(),  # type: ignore[dict-item]
        "ci_width": held["width"].median(),  # type: ignore[dict-item]
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
    parser.add_argument("--bootstrap", type=int, default=0, help="replicates for intervals")
    parser.add_argument("--configs", nargs="+", default=list(CONFIGS), choices=list(CONFIGS))
    parser.add_argument("--out", type=Path, default=Path("sim-results"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    members, units = reference(args)
    members.write_parquet(args.out / "members.parquet")
    families = units.select("unit", "family")
    rows = []
    for config in args.configs:
        build_index(args.out / "members.parquet", args.out / f"index_{config}", CONFIGS[config])
        index = Index.load(args.out / f"index_{config}")
        stats = index.meta["stats"]
        index_mb = (stats["tier1_bytes"] + stats["tier2_bytes"] + stats.get("dense_bytes", 0)) / 1e6
        for seed in range(1, args.seeds + 1):
            reads = args.out / f"reads_{seed}.fa"
            truth = sample(units, args, seed, reads)
            result = profile(index, reads, bootstrap=args.bootstrap).with_columns(
                unit=pl.col("cluster_rep").cast(pl.Int64)
            )
            for rule, (count, abundance) in RULES.items():
                found = result.filter(pl.col(count) >= 1).select("unit", estimate=abundance)
                extra = calibration(truth, result) if args.bootstrap and rule == "gather_zi" else {}
                rows.append({"config": config, "rule": rule, "seed": seed, "index_mb": index_mb,
                             **score(truth, found, families), **extra})  # fmt: skip
    scores = pl.DataFrame(rows, infer_schema_length=None)
    scores.write_csv(args.out / "scores.tsv", separator="\t")
    summary = (
        scores.group_by("config", "rule", maintain_order=True)
        .agg(pl.exclude("seed").mean())
        .with_columns(pl.selectors.float().round(3))
    )
    summary.write_csv(args.out / "summary.tsv", separator="\t")
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=200):
        print(summary)


if __name__ == "__main__":
    main()
