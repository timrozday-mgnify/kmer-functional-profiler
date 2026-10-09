"""Simulated abundance benchmark: exact per-unit truth for the phase-4 models, run locally.

Reference: ``--families`` full-length seed proteins (200-600 aa) from ``--seed-fasta``,
by default the MGnify sample (``scripts/fetch_mgnify_sample.py``; ``FL=0`` records
skipped). Each family has ``--paralogs`` units: the seed mutated to 85-95% identity
(``--paralog-identity``), so
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
- Near hits, on ``none`` rows (:func:`near_scores`): every hit unit that is not itself
  present but has a present strain in its family, scored by ``aai_naive`` against the
  highest identity of those strains to its centroid: ``near_n``, ``near_bias_<lo>`` and
  ``near_cover_<lo>`` (``aai_naive`` within 0.05) by truth bin (``NEAR_BINS``),
  ``near_spearman``, ``near_lower_bound`` (share flagged as lower bounds).

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
IDENTITIES = (1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7)  # strain to its unit's centroid
NEAR_BINS = ((0.9, 1.01), (0.8, 0.9), (0.6, 0.8))  # identity bins of near-hit AAI truth
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


SITE_GAMMA: float | None = None  # --site-gamma: shape of per-site rates; None = uniform
SITE_RATES: dict[int, np.ndarray] = {}  # site rates by sequence length


def mutate(seq: str, identity: float, rng: random.Random) -> str:
    """Substitute (1 - identity) of the sites: uniformly, or by gamma site rates
    (``--site-gamma``). Lengths never change, so one rate vector per length gives a family's
    paralogs, members and strains the same conserved sites (families of equal length share
    theirs)."""
    out, n = list(seq), round((1 - identity) * len(seq))
    rates = None
    if SITE_GAMMA:
        rates = SITE_RATES.setdefault(
            len(seq), np.random.default_rng(len(seq)).gamma(SITE_GAMMA, 1.0, len(seq))
        )
    sites = (
        rng.sample(range(len(seq)), n)
        if rates is None
        else np.random.default_rng(rng.randrange(2**32)).choice(
            len(seq), n, replace=False, p=rates / rates.sum()
        )
    )
    for i in sites:
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
            add(family, mutate(seed, rng.uniform(*args.paralog_identity), rng), None)
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
    """Write the sample's reads to ``fasta`` (FASTQ, Phred from ``--error``, if it ends in
    ``.fq``); return truth (``unit``, ``depth``, ``identity``, ``strain``, and the read bases
    in its CDS and in its flanks, ``cds_bases`` and ``flank_bases``). Strains are at
    ``--identity`` if set, else drawn from ``IDENTITIES``."""
    rng = random.Random(seed)
    identities = (args.identity,) if getattr(args, "identity", None) else IDENTITIES
    chosen = rng.sample(range(units.height), round(args.present * units.height))
    truth, reads = [], []
    for unit in chosen:
        identity = rng.choice(identities)
        depth = min(math.exp(rng.gauss(1.0, 1.0)), 50.0)
        strain = mutate(units["centroid"][unit], identity, rng)
        cds = "".join(rng.choice(SYNONYMS[aa]) for aa in strain) + "TAA"
        seq = "".join(rng.choices(BASES, k=FLANK)) + cds + "".join(rng.choices(BASES, k=FLANK))
        n = round(depth * len(seq) / READ)
        cds_bases = 0
        for _ in range(n):
            start = rng.randrange(len(seq) - READ + 1)
            cds_bases += max(0, min(start + READ, FLANK + len(cds)) - max(start, FLANK))
            read = seq[start : start + READ]
            if rng.random() < 0.5:
                read = reverse_complement(read.encode()).decode()
            reads.append(read)
        truth.append((unit, depth, identity, strain, cds_bases, n * READ - cds_bases))
    n_decoys = round(args.decoys * len(reads))
    reads += ["".join(rng.choices(BASES, k=READ)) for _ in range(n_decoys)]
    reads = [
        "".join(rng.choice(BASES.replace(b, "")) if rng.random() < args.error else b for b in r)
        for r in reads
    ]
    if fasta.suffix == ".fq":
        qual = chr(33 + round(-10 * math.log10(max(args.error, 1e-4))))
        fasta.write_text("".join(f"@r{i}\n{r}\n+\n{qual * len(r)}\n" for i, r in enumerate(reads)))
    else:
        fasta.write_text("".join(f">r{i}\n{r}\n" for i, r in enumerate(reads)))
    schema = ["unit", "depth", "identity", "strain", "cds_bases", "flank_bases"]
    return pl.DataFrame(truth, schema=schema, orient="row")


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
            f"bias_{i}": med / ref if med is not None and ref else None  # type: ignore[operator]
            for i in IDENTITIES[1:]
            for med in [tp.filter(pl.col("identity") == i)["ratio"].median()]
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


def near_scores(
    truth: pl.DataFrame, result: pl.DataFrame, units: pl.DataFrame
) -> dict[str, float | None]:
    """``aai_naive`` of near hits: hit units not present, against the best identity of their
    family's present strains to their centroid (what sylph's blanket ANI reports)."""
    family = dict(units.select("unit", "family").iter_rows())
    centroid = dict(units.select("unit", "centroid").iter_rows())
    strains: dict[int, list[str]] = {}
    for unit, strain in truth.select("unit", "strain").iter_rows():
        strains.setdefault(family[unit], []).append(strain)
    present = set(truth["unit"])
    rows = [
        (naive, low, max(_identity(st, centroid[u]) for st in strains[family[u]]))
        for u, naive, low in result.select("unit", "aai_naive", "aai_naive_lower_bound").iter_rows()
        if u not in present and family.get(u) in strains
    ]
    out: dict[str, float | None] = {"near_n": len(rows)}
    if not rows:
        return out
    df = pl.DataFrame(rows, schema=["naive", "low", "true"], orient="row")
    out["near_spearman"] = (
        df.select(pl.corr("naive", "true", method="spearman")).item() if df.height > 2 else None
    )
    out["near_lower_bound"] = df["low"].mean()  # type: ignore[assignment]
    for lo, hi in NEAR_BINS:
        part = df.filter(pl.col("true").is_between(lo, hi, closed="left"))
        out[f"near_bias_{lo}"] = (part["naive"] - part["true"]).median() if part.height else None  # type: ignore[assignment]
        out[f"near_cover_{lo}"] = (
            ((part["naive"] - part["true"]).abs() <= 0.05).mean() if part.height else None
        )
    return out


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


def unknown_ladder(args: argparse.Namespace, members: pl.DataFrame, units: pl.DataFrame) -> None:
    """Phase 7, step 14's evaluation of the unknown fraction: per ``--holdout`` fraction h,
    the first h of the families (in a fixed shuffled order, so the ladder is nested) are left
    out of the index, whole, so their reads have no relative in it. Each seed's sample (FASTQ,
    so ``error_thinning`` sees ``--error``) is profiled with ``summary``; truth is the share of
    read bases in CDSs of indexed units (``truth_known``), of held-out units
    (``truth_heldout``), in flanks (``truth_noncoding``) and decoys (``truth_decoy``).
    ``census_containment`` is read against them (``explained_fraction`` was dropped in
    step 40). Writes ``unknown_scores.tsv`` and ``unknown_summary.tsv``."""
    order = list(range(int(units["family"].max()) + 1))  # type: ignore[arg-type]
    random.Random(1).shuffle(order)
    rows = []
    for config, h in itertools.product(args.configs, args.holdout):
        held = set(order[: round(h * len(order))])
        kept = units.filter(~pl.col("family").is_in(held))["unit"]
        name = f"{config}_h{round(h * 100)}"
        members.filter(pl.col("cluster_rep").is_in(kept.to_list())).write_parquet(
            args.out / "kept.parquet"
        )
        build_index(args.out / "kept.parquet", args.out / f"index_{name}", CONFIGS[config])
        index = Index.load(args.out / f"index_{name}")
        for seed in range(1, args.seeds + 1):
            reads = args.out / f"reads_{seed}.fq"
            truth = sample(units, args, seed, reads).join(units.select("unit", "family"), on="unit")
            summary: dict[str, float | int | None] = {}
            profile(index, reads, summary=summary)
            bases = summary["bases"]
            assert bases
            is_held = pl.col("family").is_in(held)
            cds = truth.select(
                known=pl.col("cds_bases").filter(~is_held).sum(),
                heldout=pl.col("cds_bases").filter(is_held).sum(),
                flanks=pl.col("flank_bases").sum(),
            ).row(0, named=True)
            row = {
                "config": config, "holdout": h, "identity": args.identity, "seed": seed,
                "truth_known": cds["known"] / bases, "truth_heldout": cds["heldout"] / bases,
                "truth_noncoding": cds["flanks"] / bases,
                "truth_decoy": 1 - (cds["known"] + cds["heldout"] + cds["flanks"]) / bases,
                **{k: summary[k] for k in ("census_containment", "error_thinning",
                                           "census_kmers")},
            }  # fmt: skip
            rows.append(row)
    scores = pl.DataFrame(rows, infer_schema_length=None)
    scores.write_csv(args.out / "unknown_scores.tsv", separator="\t")
    summary_df = (
        scores.group_by("config", "holdout", maintain_order=True)
        .agg(pl.exclude("seed", "identity").mean())
        .with_columns(pl.selectors.float().round(3))
    )
    summary_df.write_csv(args.out / "unknown_summary.tsv", separator="\t")
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=200):
        print(summary_df)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-fasta", type=Path, default=SEEDS_FASTA)
    parser.add_argument("--families", type=int, default=100)
    parser.add_argument("--paralogs", type=int, default=3)
    parser.add_argument("--members", type=int, default=4)
    parser.add_argument(
        "--paralog-identity",
        type=float,
        nargs=2,
        default=(0.85, 0.95),
        metavar=("LO", "HI"),
        help="paralog centroids' identity to the seed",
    )
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
    parser.add_argument(
        "--holdout",
        type=float,
        nargs="+",
        help="run the unknown-fraction ladder: shares of families held out",
    )
    parser.add_argument("--identity", type=float, help="every strain at this identity")
    parser.add_argument(
        "--site-gamma",
        type=float,
        help="gamma shape of per-site substitution rates (clustered variation)",
    )
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    global SITE_GAMMA
    SITE_GAMMA = args.site_gamma

    members, units = reference(args)
    members.write_parquet(args.out / "members.parquet")
    if args.holdout:
        unknown_ladder(args, members, units)
        return
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
                    else near_scores(truth, result, units)
                    if rule == "none"
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
