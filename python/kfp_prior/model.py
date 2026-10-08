"""Genome-informed unit presence (plan: Genome-informed unit presence, option B1).

A companion of the profiler, as Bracken is to Kraken: it reads the profiler's and genome
mode's output files (``profile.tsv``, ``genomes.tsv``, the genome index of
``annotate-genomes``) and never imports the profiler.

- :class:`Carriage`: carriage counts per family, genus and species and a shrinkage strength
  α per rank, as the profiler's ``species-index --genomes`` writes them (the counts and α
  fit moved there with the species model, phase 11, step 7).
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
Q_FLOOR: Final = 0.05  # carriage frequencies below this give no prior
CLIP: Final = 1e-6


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Carriage:
    """The carriage table of a species index built from a genome index
    (``species-index --genomes``): ``carriage.parquet``, ``clades.parquet``,
    ``lineage.parquet`` and α in ``meta.json``."""

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
    if _sha256(genome_index / "meta.json") != carriage.meta.get("genome_index_meta_sha256"):
        raise ValueError(
            f"the carriage table was built from another genome index than {genome_index}"
        )
    n_units = json.loads((genome_index / "meta.json").read_text())["units"]
    observed = profile.select(
        pl.col("unit").cast(pl.UInt32), pl.col("hits").cast(pl.Float64), "present_prob"
    ).filter(pl.col("hits") > 0)
    expected_present = float(observed["present_prob"].sum())
    # floored at one unit: a profile whose present_prob sums to 0 keeps its zeros (0 odds)
    pi0 = min(max(expected_present, 1.0) / n_units, 1 - CLIP)
    numeric = [c for c in ("depth", "present", "present_prob") if c in genomes.columns]
    weight = (
        genomes.cast(dict.fromkeys(numeric, pl.Float64))
        .select(
            pl.col("genome").cast(pl.UInt32),
            rate=pl.col("depth") * pl.col("present"),
            p_g=pl.col("present_prob") if "present_prob" in genomes.columns else pl.lit(1.0),
        )
        .filter(pl.col("rate") > 0)
    )
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
