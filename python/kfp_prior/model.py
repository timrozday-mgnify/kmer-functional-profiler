"""Genome-informed unit presence (plan: Genome-informed unit presence, option B1).

A companion of the profiler, as Bracken is to Kraken: it reads the profiler's and genome
mode's output files (``profile.tsv``, ``genomes.tsv``, the genome index of
``annotate-genomes``) and never imports the profiler.

- :func:`build_carriage`: from which units each reference genome carries (its proteins' best
  units, ``genome_best.parquet``), carriage counts per family, genus and species, and a
  shrinkage strength α per rank fitted by leave-one-genome-out log loss.
- :func:`update`: per detected genome G, the carriage frequency q_{G,u} of its species
  shrunk towards genus, family and all genomes; a noisy-OR prior over genomes; and each
  unit's presence updated by its own hits (closed form at zero hits, odds rescaling of the
  profile's ``present_prob`` otherwise).
"""

import hashlib
import json
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl

RANKS: Final = ("family", "genus", "species")  # top down
GTDB_PREFIX: Final = {"f": "family", "g": "genus", "s": "species"}
ALPHAS: Final = np.geomspace(0.01, 100, 21)
MAX_HOLDOUT: Final = 200  # held-out genomes per rank when fitting α
Q_FLOOR: Final = 0.05  # carriage frequencies below this give no prior
CLIP: Final = 1e-6


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lineages(genomes: pl.DataFrame) -> pl.DataFrame:
    """``genome`` and its ``family``, ``genus`` and ``species`` from a GTDB-style
    ``taxonomy`` column (by prefix; null where missing). A genome without a species is a
    species of its own (``genome:<name>``)."""
    taxa = genomes["taxonomy"] if "taxonomy" in genomes.columns else [None] * genomes.height
    rows = []
    for g, name, taxonomy in zip(genomes["genome"], genomes["name"], taxa, strict=True):
        parts = [x.strip() for x in (taxonomy or "").split(";")]
        by_rank = {GTDB_PREFIX[x[0]]: x for x in parts if x[1:3] == "__" and len(x) > 3
                   and x[0] in GTDB_PREFIX}  # fmt: skip
        rows.append((g, by_rank.get("family"), by_rank.get("genus"),
                     by_rank.get("species") or f"genome:{name}"))  # fmt: skip
    schema = {"genome": pl.UInt32, **dict.fromkeys(RANKS, pl.String)}
    return pl.DataFrame(rows, schema=schema, orient="row")


def _counts(carried: pl.DataFrame, lineage: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Carriage counts ``rank``, ``clade``, ``unit``, ``n`` and clade sizes ``rank``,
    ``clade``, ``N``, per rank and for all genomes (rank ``root``, clade "")."""
    long = pl.concat(
        [lineage.select("genome", rank=pl.lit("root"), clade=pl.lit(""))]
        + [lineage.select("genome", rank=pl.lit(r), clade=pl.col(r)) for r in RANKS]
    ).drop_nulls("clade")
    n = carried.join(long, on="genome").group_by("rank", "clade", "unit").agg(n=pl.len())
    sizes = long.group_by("rank", "clade").agg(N=pl.len())
    return n.cast({"n": pl.Int64}), sizes.cast({"N": pl.Int64})


def _parent(lineage: pl.DataFrame, rank: str) -> pl.DataFrame:
    """Per genome, the nearest clade above ``rank`` (``p_rank``, ``p``), else root."""
    above = RANKS[: RANKS.index(rank)]
    p_rank, p = pl.lit("root"), pl.lit("")
    for r in above:  # top down: the last non-null wins
        p_rank = pl.when(pl.col(r).is_not_null()).then(pl.lit(r)).otherwise(p_rank)
        p = pl.when(pl.col(r).is_not_null()).then(pl.col(r)).otherwise(p)
    return lineage.select("genome", p_rank=p_rank, p=p)


def fit_alpha(
    carried: pl.DataFrame, lineage: pl.DataFrame, rank: str, seed: int = 0
) -> tuple[float, dict[str, float]]:
    """α at ``rank`` by leave-one-genome-out log loss: each held-out genome's carriage of
    the units its parent clade carries (its own clade's, when the parent is root), predicted
    from the rest of its clade shrunk towards the rest of the parent's frequency. Returns
    α and the losses of the best α, of no shrinkage (α → 0) and of the parent alone
    (α → ∞), per prediction."""
    n, sizes = _counts(carried, lineage)
    eligible = (
        lineage.select("genome", clade=pl.col(rank))
        .join(sizes.filter(pl.col("rank") == rank), on="clade")
        .filter(pl.col("N") >= 2)
    )
    if eligible.height == 0:
        return 1.0, {}
    held = eligible.sample(min(MAX_HOLDOUT, eligible.height), seed=seed).join(
        _parent(lineage, rank), on="genome"
    )
    by_parent = held.filter(pl.col("p_rank") != "root").join(
        n, left_on=["p_rank", "p"], right_on=["rank", "clade"]
    )
    by_clade = held.filter(pl.col("p_rank") == "root").join(
        n.filter(pl.col("rank") == rank), on="clade"
    )
    rows = (
        pl.concat(
            [
                by_parent.select("genome", "clade", "p_rank", "p", "unit"),
                by_clade.select("genome", "clade", "p_rank", "p", "unit"),
            ]
        )  # fmt: skip
        .join(
            n.filter(pl.col("rank") == rank).select("clade", "unit", n_c="n"),
            on=["clade", "unit"],
            how="left",
        )  # fmt: skip
        .join(sizes.filter(pl.col("rank") == rank).select("clade", N_c="N"), on="clade")
        .join(
            n.select(p_rank="rank", p="clade", unit="unit", n_p="n"),
            on=["p_rank", "p", "unit"],
            how="left",
        )  # fmt: skip
        .join(sizes.select(p_rank="rank", p="clade", N_p="N"), on=["p_rank", "p"])
        .join(carried.with_columns(x=pl.lit(1)), on=["genome", "unit"], how="left")
        .with_columns(pl.col("n_c", "n_p", "x").fill_null(0))
    )
    x = rows["x"].to_numpy().astype(np.float64)
    rest = rows["n_c"].to_numpy() - x
    parent = (rows["n_p"].to_numpy() - x) / np.maximum(rows["N_p"].to_numpy() - 1, 1)
    size = rows["N_c"].to_numpy() - 1.0

    def loss(alpha: float) -> float:
        q = np.clip((rest + alpha * parent) / (size + alpha), CLIP, 1 - CLIP)
        return float(-(x * np.log(q) + (1 - x) * np.log1p(-q)).mean())

    losses = [loss(a) for a in ALPHAS]
    best = float(ALPHAS[int(np.argmin(losses))])
    return best, {"best": min(losses), "no_shrinkage": loss(1e-9), "parent_only": loss(1e9)}


def build_carriage(genome_index: str | Path, out: str | Path) -> dict[str, object]:
    """Write the carriage table of a genome index (``annotate-genomes`` output, with
    ``genome_best.parquet`` and a ``taxonomy`` column in ``genomes.tsv`` for shrinkage):
    ``carriage.parquet`` (``rank``, ``clade``, ``unit``, ``n``), ``clades.parquet``
    (``rank``, ``clade``, ``N``), ``lineage.parquet`` and ``meta.json`` (α per rank, its
    held-out losses, the genome index's checksum)."""
    genome_index, out = Path(genome_index), Path(out)
    genomes = pl.read_csv(genome_index / "genomes.tsv", separator="\t")
    lineage = lineages(genomes.cast({"genome": pl.UInt32}))
    carried = (
        pl.read_parquet(genome_index / "genome_best.parquet").select("genome", "unit").unique()
    )
    n, sizes = _counts(carried, lineage)
    fits = {r: fit_alpha(carried, lineage, r) for r in RANKS}
    out.mkdir(parents=True, exist_ok=True)
    n.sort("rank", "clade", "unit").write_parquet(out / "carriage.parquet")
    sizes.sort("rank", "clade").write_parquet(out / "clades.parquet")
    lineage.write_parquet(out / "lineage.parquet")
    meta = {
        "genome_index_meta_sha256": _sha256(genome_index / "meta.json"),
        "alpha": {r: a for r, (a, _) in fits.items()},
        "held_out_log_loss": {r: losses for r, (_, losses) in fits.items()},
        "genomes": genomes.height,
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    return meta


class Carriage:
    """A carriage table written by :func:`build_carriage`."""

    def __init__(self, path: str | Path) -> None:
        path = Path(path)
        self.meta = json.loads((path / "meta.json").read_text())
        self.n = pl.read_parquet(path / "carriage.parquet")
        self.sizes = pl.read_parquet(path / "clades.parquet")
        self.lineage = pl.read_parquet(path / "lineage.parquet")

    def q(self, genomes: list[int]) -> pl.DataFrame:
        """Per detected genome and candidate unit, ``q``: the frequency of carrying the unit
        in its species, shrunk towards genus, family and all genomes (α per rank), for the
        units its highest known clade carries; rows with q < ``Q_FLOOR`` are dropped."""
        parts = []
        root = self.n.filter(pl.col("rank") == "root").select("unit", n_root="n")
        n_root = self.sizes.filter(pl.col("rank") == "root")["N"][0]
        for g in genomes:
            row = self.lineage.filter(pl.col("genome") == g).row(0, named=True)
            chain = [(r, row[r]) for r in RANKS if row[r] is not None]
            top_rank, top = chain[0]
            units = self.n.filter((pl.col("rank") == top_rank) & (pl.col("clade") == top))
            q = (
                units.select("unit")
                .join(root, on="unit")
                .select("unit", q=pl.col("n_root") / n_root)
            )
            for r, clade in chain:
                counts = self.n.filter((pl.col("rank") == r) & (pl.col("clade") == clade))
                size = self.sizes.filter((pl.col("rank") == r) & (pl.col("clade") == clade))
                alpha = self.meta["alpha"][r]
                q = q.join(counts.select("unit", "n"), on="unit", how="left").select(
                    "unit",
                    q=(pl.col("n").fill_null(0) + alpha * pl.col("q")) / (size["N"][0] + alpha),
                )
            parts.append(q.filter(pl.col("q") >= Q_FLOOR).with_columns(genome=pl.lit(g)))
        schema = pl.Schema({"unit": pl.UInt32, "q": pl.Float64, "genome": pl.UInt32})
        return pl.concat([pl.DataFrame(schema=schema), *(p.cast(schema) for p in parts)])


def zero_hit_presence(rho: np.ndarray, rate: np.ndarray) -> np.ndarray:
    """P(carriage via a genome | 0 hits) for prior ``rho`` and expected hits ``rate``."""
    kept = rho * np.exp(-rate)
    out: np.ndarray = kept / (1 - rho + kept)
    return out


def _content(genome_index: Path, units: np.ndarray) -> pl.DataFrame:
    """``unit``, ``genome``, ``c`` rows of ``units`` from the genome index's unit-major
    CSR (the file contract of ``annotate-genomes``)."""
    offsets = np.load(genome_index / "unit_offsets.npy", mmap_mode="r")
    genome = np.load(genome_index / "genome.npy", mmap_mode="r")
    hits = np.load(genome_index / "hits.npy", mmap_mode="r")
    units = np.asarray(units, dtype=np.int64)
    starts = offsets[units].astype(np.int64)
    lengths = offsets[units + 1].astype(np.int64) - starts
    rows = np.concatenate([np.arange(s, s + n) for s, n in zip(starts, lengths, strict=True)]
                          + [np.empty(0, np.int64)])  # fmt: skip
    return pl.DataFrame(
        {
            "unit": np.repeat(units, lengths).astype(np.uint32),
            "genome": np.asarray(genome[rows], dtype=np.uint32),
            "c": np.asarray(hits[rows], dtype=np.float64),
        }  # fmt: skip
    )


def _odds(p: pl.Expr) -> pl.Expr:
    return p / (1 - p)


def update(
    profile: pl.DataFrame,
    genomes: pl.DataFrame,
    genome_index: str | Path,
    carriage: Carriage,
) -> tuple[pl.DataFrame, pl.DataFrame | None]:
    """Unit presence updated with a genome-informed prior.

    ``genomes`` is genome mode's ``genomes.tsv`` (``genome``, ``depth`` λ_G, ``present``
    π_G; ``present_prob`` P_G if given, else 1). Each detected genome carries unit *u* with
    prior ρ_{G,u} = P_G q_{G,u} (:meth:`Carriage.q`), giving expected hits
    λ_G π_G e_{G,u} (e: G's own content on *u* if it carries it, else the mean content of
    the genomes that do). The prior is a noisy-OR over genomes and the profile's global
    prior π₀ = Σ ``present_prob`` / index units (the odds ``present_prob`` was fitted at):
    1 - (1 - π₀) Π_G (1 - ρ_{G,u}).

    - Hit units (in the profile): ``present_prob``'s odds x prior odds / π₀'s odds, so a
      unit no detected genome carries keeps ``present_prob`` exactly.
    - Units with zero hits: genomes factorise, P = 1 - Π_G (1 - ρ e^(-r) / (1 - ρ + ρ e^(-r)))
      with r the expected hits; ``imputed`` if ≥ 0.5. ``sibling_hit``: one of its Pfams is
      observed through another unit (the function is there; this unit may not be).

    Returns ``presence`` (``unit``, ``hits``, ``present_prob``, ``prior``, ``expected_hits``,
    ``present_prob_updated``, ``imputed``, ``sibling_hit``) and, if the genome index has
    Pfam labels, per Pfam ``present_prob_observed`` and ``present_prob_updated`` (units
    taken as independent) and ``imputed``.
    """
    genome_index = Path(genome_index)
    if _sha256(genome_index / "meta.json") != carriage.meta["genome_index_meta_sha256"]:
        raise ValueError(
            f"the carriage table was built from another genome index than {genome_index}"
        )
    n_units = json.loads((genome_index / "meta.json").read_text())["units"]
    observed = profile.select(
        pl.col("unit").cast(pl.UInt32), pl.col("hits").cast(pl.Float64), "present_prob"
    ).filter(pl.col("hits") > 0)
    expected_present = float(observed["present_prob"].sum())
    pi0 = expected_present / n_units
    weight = genomes.select(
        pl.col("genome").cast(pl.UInt32),
        rate=pl.col("depth") * pl.col("present"),
        p_g=pl.col("present_prob") if "present_prob" in genomes.columns else pl.lit(1.0),
    ).filter(pl.col("rate") > 0)
    q = carriage.q(weight["genome"].to_list())
    content = _content(genome_index, np.unique(q["unit"].to_numpy()))
    best = pl.read_parquet(genome_index / "genome_best.parquet").select("genome", "unit").unique()
    carried = content.join(best, on=["genome", "unit"], how="semi")
    typical = carried.group_by("unit").agg(e_mean=pl.col("c").mean())
    pairs = (
        q.join(weight, on="genome")
        .join(carried.rename({"c": "e_own"}), on=["unit", "genome"], how="left")
        .join(typical, on="unit", how="left")
        .with_columns(
            rho=pl.col("p_g") * pl.col("q"),
            hits_expected=pl.col("rate") * pl.coalesce("e_own", "e_mean").fill_null(0.0),
        )
    )
    rho, rate = pairs["rho"].to_numpy(), pairs["hits_expected"].to_numpy()
    pairs = pairs.with_columns(
        log_not=np.log1p(-np.minimum(rho, 1 - CLIP)),
        log_not_zero=np.log1p(-np.minimum(zero_hit_presence(rho, rate), 1 - CLIP)),
        weighted=rho * rate,
    )
    per_unit = pairs.group_by("unit").agg(
        prior=1 - (1 - pi0) * pl.col("log_not").sum().exp(),
        expected_hits=pl.col("weighted").sum(),
        zero=1 - pl.col("log_not_zero").sum().exp(),
    )
    prior_odds = _odds(pl.col("prior").fill_null(pi0)) / (pi0 / (1 - pi0))
    hit = (
        observed.join(per_unit, on="unit", how="left")
        .with_columns(
            pl.col("prior").fill_null(pi0),
            pl.col("expected_hits").fill_null(0.0),
            present_prob_updated=pl.when(pl.col("present_prob") >= 1)
            .then(1.0)
            .otherwise(
                (_odds(pl.col("present_prob")) * prior_odds)
                / (1 + _odds(pl.col("present_prob")) * prior_odds)
            ),
        )
        .drop("zero")
    )
    zero = (
        per_unit.join(observed, on="unit", how="anti").select(
            "unit",
            hits=pl.lit(0.0),
            present_prob=pl.lit(None, pl.Float64),
            prior="prior",
            expected_hits="expected_hits",
            present_prob_updated="zero",
        )  # fmt: skip
    )
    presence = pl.concat([hit, zero], how="diagonal_relaxed").with_columns(
        imputed=(pl.col("hits") == 0) & (pl.col("present_prob_updated") >= 0.5)
    )
    pfam_path = genome_index / "unit_pfam.parquet"
    if not pfam_path.exists():
        return presence.with_columns(sibling_hit=pl.lit(None, pl.Boolean)).sort("unit"), None
    labels = pl.read_parquet(pfam_path).select(pl.col("unit").cast(pl.UInt32), "pfam_accession")
    seen = (
        labels.join(
            presence.filter((pl.col("hits") > 0) & (pl.col("present_prob") >= 0.5)), on="unit"
        )["pfam_accession"]
        .unique()
        .implode()
    )
    sibling = labels.filter(pl.col("pfam_accession").is_in(seen)).select("unit").unique()
    presence = presence.with_columns(
        sibling_hit=(pl.col("hits") == 0) & pl.col("unit").is_in(sibling["unit"].implode())
    ).sort("unit")
    by_pfam = (
        presence.join(labels, on="unit")
        .group_by("pfam_accession")
        .agg(
            present_prob_observed=1 - (1 - pl.col("present_prob").fill_null(0.0)).product(),
            present_prob_updated=1 - (1 - pl.col("present_prob_updated")).product(),
        )
        .with_columns(
            imputed=(pl.col("present_prob_observed") < 0.5)
            & (pl.col("present_prob_updated") >= 0.5)
        )
        .sort("pfam_accession")
    )
    return presence, by_pfam
