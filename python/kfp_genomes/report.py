"""``place``'s outputs (phase 12, step 5): the genome profile, posterior samples, per-strain
unit carriage, model checks and model-free evidence.

Per detected species (P(K ≥ 1) ≥ ``DETECTED``) the number of strains is the posterior mode
over K = 1 .. ``max_strains`` (exact for K = 1, stepping-stone above); its strains are
summarised from posterior draws:

- depth (mean, 95% interval) and relative abundance (over every reported strain);
- placement: the commonest (g, h, t, ℓ) of the draws; the nearest references with posterior
  weights (a draw between g and h at t gives g 1 - t and h t); ℓ's mean; ``novel`` when
  P(ℓ > 0.5) > 0.5;
- species ``present_prob`` = P(K ≥ 1) and ``mixture_prob`` = P(K ≥ 2).

*Model checks* (posterior predictive, reported, never patched; a species fails one at
p < ``CHECK_ALPHA``):

- ``check_core``: core units (q ≥ 0.9) detected against the number predicted, P(carried |
  placement) x P(≥ 1 hit | depth, survival) per unit (Poisson-binomial, normal
  approximation);
- ``check_spread``: the depths of hit units the strains carry against predicted at the
  references' copies (z on log coverage, its sd from the likelihood-ratio interval
  ``coverage_interval`` gives, rebuilt from the evidence as the profile's own): units with
  |z| > ``OUTLIER_Z`` against the binomial rate copy changes (``COPY_CHANGE``) and chance
  allow. (A χ² over all units, as first written, failed every species with a few copy
  changes.)
- ``check_single``: for K ≥ 2, the depth of units only one strain carries against that
  strain's (mean z per strain, the smallest p).

*Model-free evidence* (no fit; to read beside the model's numbers): over the hit units
(``present_prob`` ≥ 0.5) of the detected species, ``units_hit`` and ``core_hit`` of
``core_units``; ``units_own``, the units gather gives the species when each hit unit goes to
the first species, in gather order, that has it (G0, one level up); ``kmers_own``, those
units' ``kmers_unique``; and ``depth_greedy``, the median ``coverage_zi`` / copies over its
own core units (all its own units if fewer than ``MIN_GREEDY_UNITS``), sylph-like.
"""

from dataclasses import dataclass
from typing import Final

import numpy as np
import polars as pl
from scipy.stats import binom, norm

from kfp_genomes.evidence import Evidence, Histogram, log_likelihood
from kfp_genomes.mixtures import fit_mixture, strain_number_posterior
from kfp_genomes.panel import Panel, SpeciesPanel
from kfp_genomes.place import (
    COPY_CHANGE,
    CORE,
    NOVEL,
    OneStrain,
    Placed,
    place_species,
    survival,
)
from kmer_functional_profiler.query import coverage_interval, gather, ztp_lambda

DETECTED: Final = 0.5
CHECK_ALPHA: Final = 0.01
NEAREST: Final = 3
DRAWS: Final = 500  # posterior draws kept per K = 1 strain
CHECK_DRAWS: Final = 100  # draws the predictive checks average over
MIN_GREEDY_UNITS: Final = 5
OUTLIER_Z: Final = 3.0  # |z| of a unit's log coverage from which it is an outlier


@dataclass
class Strains:
    """One species' reported strains: draws (draws x K points and λ), per-strain unit
    carriage and depth (K x U), the K posterior and diagnostics."""

    points: np.ndarray
    lam: np.ndarray
    carriage: np.ndarray
    depth: np.ndarray
    survival: np.ndarray  # U: the strains' survival (the species' mean x scale)
    k_post: dict[int, float]
    diagnostics: dict[str, float]


def _strains(
    placed: Placed,
    s: int,
    max_strains: int,
    hist: Histogram | None,
    rng: np.random.Generator,
    options: dict[str, int],
) -> Strains:
    one, sp, terms = placed.fits[s], placed.panels[s], placed.terms[s]
    log_bf = {1: one.log_bf}
    mixtures = {}
    for k in range(2, max_strains + 1):
        mixtures[k] = fit_mixture(sp, terms, one, k, hist=hist, seed=int(rng.integers(2**31)),
                                  **options)  # type: ignore[arg-type]  # fmt: skip
        log_bf[k] = mixtures[k].log_bf
    k_post = strain_number_posterior(log_bf, placed.pi)
    k_map = max((k for k in k_post if k >= 1), key=lambda k: k_post[k])
    scale = float(one.extra["survival_scale"])
    surv = survival(sp, scale)[-1]
    if k_map == 1:
        point, lam = one.sample(DRAWS, rng)
        return Strains(point[:, None], lam[:, None], one.unit_carriage[None],
                       one.unit_depth[None], surv, k_post, {})  # fmt: skip
    mix = mixtures[k_map]
    scale = mix.diagnostics.get("survival_scale", scale)
    return Strains(mix.points, mix.lam, mix.carriage, mix.depth, survival(sp, scale)[-1],
                   k_post, mix.diagnostics)  # fmt: skip


def _placement(sp: SpeciesPanel, one: OneStrain, points: np.ndarray) -> dict[str, object]:
    """Summaries of one strain's draws of placement points."""
    grid = one.grid
    g, h, t, ell = grid.g[points], grid.h[points], grid.t[points], grid.ell[points]
    weight = np.bincount(g, (1 - t) * (1 - ell), len(sp.genomes)) + np.bincount(
        h, t * (1 - ell), len(sp.genomes)
    )
    top = np.argsort(-weight)[:NEAREST]
    total = weight.sum()
    values, counts = np.unique(points, return_counts=True)
    mode = int(values[np.argmax(counts)])
    return {
        "ref_g": sp.names[grid.g[mode]],
        "ref_h": sp.names[grid.h[mode]],
        "t": float(grid.t[mode]),
        "ell": float(grid.ell[mode]),
        "ell_mean": float(ell.mean()),
        "novel": bool((ell > NOVEL).mean() > 0.5),
        "nearest": ";".join(
            f"{sp.names[i]}:{weight[i] / total:.3f}" for i in top if total > 0 and weight[i] > 0
        ),  # fmt: skip
    }


def _carry_prior(sp: SpeciesPanel, one: OneStrain, points: np.ndarray) -> np.ndarray:
    """P(carried | placement) per draw and unit: draws x U."""
    grid = one.grid
    g, h, t, ell = grid.g[points], grid.h[points], grid.t[points], grid.ell[points]
    return np.asarray(
        (1 - ell)[:, None] * ((1 - t)[:, None] * sp.xt[g] + t[:, None] * sp.xt[h])
        + ell[:, None] * sp.q[None]
    )


def model_checks(placed: Placed, s: int, st: Strains) -> dict[str, float]:
    """The posterior predictive checks of the module docstring: p-values."""
    one, sp, terms = placed.fits[s], placed.panels[s], placed.terms[s]
    out: dict[str, float] = {}
    pick = np.linspace(0, len(st.lam) - 1, min(CHECK_DRAWS, len(st.lam))).astype(int)
    core = sp.q >= CORE
    if core.any():
        miss = np.ones((len(pick), int(core.sum())))
        for k in range(st.lam.shape[1]):
            carry = _carry_prior(sp, one, st.points[pick, k])[:, core]
            depth = st.lam[pick, k][:, None] * sp.n_mean[core][None]
            p_hit = -np.expm1(log_likelihood(0, 0, sp.m_g[core][None], depth,
                                             st.survival[core][None]))  # fmt: skip
            miss *= 1 - carry * p_hit
        p = 1 - miss.mean(0)
        observed = float((terms.hit | terms.blank)[core].sum())
        z = (observed - p.sum()) / np.sqrt(max(float((p * (1 - p)).sum()), 1e-9))
        out["check_core"] = float(2 * norm.sf(abs(z)))
    carried = st.carriage.sum(0)
    # predicted depth at the references' copies (a copy change is an outlier here)
    copies = np.average(np.vstack([sp.n, sp.n_mean[None]]), axis=0)
    expect = np.zeros(len(sp.units))
    for k in range(st.lam.shape[1]):
        expect += st.carriage[k] * float(np.mean(st.lam[:, k])) * copies
    use = terms.hit & (carried >= 0.9) & (expect > 0)
    if use.any():
        # coverage_zi and its likelihood-ratio interval, rebuilt as the profile has them
        c_hat = np.maximum(ztp_lambda(terms.h[use] / np.maximum(terms.k[use], 1e-9)), 1e-9)
        lo, hi = coverage_interval(c_hat, terms.k[use], terms.w[use])
        z_975 = float(norm.ppf(0.975))
        sd = np.where(lo > 0, np.log(hi / np.maximum(lo, 1e-12)) / (2 * z_975),
                      np.log(hi / c_hat) / z_975)  # fmt: skip
        z = (np.log(c_hat) - np.log(expect[use])) / np.maximum(sd, 1e-6)
        outliers = int((np.abs(z) > OUTLIER_Z).sum())
        rate = COPY_CHANGE + 2 * norm.sf(OUTLIER_Z)
        out["check_spread"] = float(binom.sf(outliers - 1, int(use.sum()), rate))
        if st.lam.shape[1] > 1:
            ps = []
            for k in range(st.lam.shape[1]):
                alone = (st.carriage[k] >= 0.9) & (np.delete(st.carriage, k, 0).max(0) <= 0.1)
                sel = alone[use]
                if sel.sum() >= 3:
                    zk = z[sel]
                    ps.append(float(2 * norm.sf(abs(zk.mean()) * np.sqrt(len(zk)))))
            if ps:
                out["check_single"] = min(ps)
    return out


def greedy_evidence(
    profile: pl.DataFrame, ev: Evidence, panel: Panel, species: list[int]
) -> pl.DataFrame:
    """The model-free columns of the module docstring, per species in ``species``."""
    hit = pl.DataFrame({"unit": ev.unit[ev.present_prob >= 0.5].astype(np.uint32)})
    su = panel.species_units.filter(pl.col("species").is_in(species))
    empty = pl.DataFrame(
        {"species": species},
        schema={"species": pl.UInt32},
    ).with_columns(
        units_hit=pl.lit(0), core_units=pl.lit(0), core_hit=pl.lit(0), units_own=pl.lit(0),
        kmers_own=pl.lit(0.0), depth_greedy=pl.lit(None, pl.Float64),
    )  # fmt: skip
    if not species:
        return empty
    hits = su.join(hit, on="unit", how="semi")
    counts = (
        su.group_by("species")
        .agg(core_units=(pl.col("q") >= CORE).sum())
        .join(
            hits.group_by("species").agg(units_hit=pl.len(), core_hit=(pl.col("q") >= CORE).sum()),
            on="species",
            how="left",
        )
    )
    own = (
        gather(
            hits.select(unit="species", hash=pl.col("unit").cast(pl.UInt64)),
            np.ones(int(max(species)) + 1),
        )  # fmt: skip
        if hits.height
        else pl.DataFrame(schema={"unit": pl.UInt32, "gather_rank": pl.UInt32})
    )
    # each hit unit's species: the first in gather order among those with it
    first = (
        hits.join(own.select(species="unit", rank="gather_rank"), on="species")
        .sort("rank")
        .unique("unit", keep="first")
    )
    cols = ["unit", "coverage_zi", *(["kmers_unique"] if "kmers_unique" in profile.columns
                                     else [])]  # fmt: skip
    per_unit = first.join(profile.select(cols).cast({"unit": pl.UInt32}), on="unit", how="left")
    if "kmers_unique" not in per_unit.columns:  # simulated profiles: the fit's hit k-mers
        k = pl.DataFrame({"unit": ev.unit.astype(np.uint32), "kmers_unique": ev.k})
        per_unit = per_unit.join(k, on="unit", how="left")
    per_unit = per_unit.join(
        su.select("species", "unit", "n_mean"), on=["species", "unit"], how="left"
    ).with_columns(per_copy=pl.col("coverage_zi") / pl.col("n_mean"))
    greedy = per_unit.group_by("species").agg(
        units_own=pl.len(),
        kmers_own=pl.col("kmers_unique").cast(pl.Float64).sum(),
        depth_core=pl.col("per_copy").filter(pl.col("q") >= CORE).median(),
        n_core=(pl.col("q") >= CORE).sum(),
        depth_all=pl.col("per_copy").median(),
    )
    return (
        counts.join(greedy, on="species", how="left")
        .select(
            pl.col("species").cast(pl.UInt32),
            pl.col("units_hit").fill_null(0),
            "core_units",
            pl.col("core_hit").fill_null(0),
            pl.col("units_own").fill_null(0),
            pl.col("kmers_own").fill_null(0.0),
            depth_greedy=pl.when(pl.col("n_core") >= MIN_GREEDY_UNITS)
            .then(pl.col("depth_core"))
            .otherwise(pl.col("depth_all")),
        )
        .sort("species")
    )


def genome_profile(
    profile: pl.DataFrame,
    panel: Panel,
    *,
    max_strains: int = 1,
    hist: Histogram | None = None,
    species: list[int] | None = None,
    seed: int = 0,
    sampler: dict[str, int] | None = None,
    **place_options: object,
) -> dict[str, object]:
    """``place``: the genome profile and its companions (module docstring). Returns
    ``genome_profile``, ``candidates``, ``placements``, ``strain_units`` (data frames) and
    ``summary``."""
    rng = np.random.default_rng(seed)
    placed = place_species(profile, panel, species=species, **place_options)  # type: ignore[arg-type]
    names = panel.species.select(pl.col("species").cast(pl.UInt32), "id", "name", "taxonomy")
    detected = sorted(s for s, f in placed.fits.items() if f.present_prob >= DETECTED)
    rows, draws, units = [], [], []
    for s in detected:
        one, sp = placed.fits[s], placed.panels[s]
        st = _strains(placed, s, max_strains, hist, rng, sampler or {})
        checks = model_checks(placed, s, st)
        failed = ",".join(k.removeprefix("check_") for k, p in checks.items() if p < CHECK_ALPHA)
        present = 1 - st.k_post.get(0, 0.0)
        mixture = sum(p for k, p in st.k_post.items() if k >= 2)
        k = st.lam.shape[1]
        for j in range(k):
            lam = st.lam[:, j]
            if k == 1:
                depth, lo, hi = one.depth_summary()
            else:
                lo, hi = (float(x) for x in np.percentile(lam, [2.5, 97.5]))
                depth = float(lam.mean())
            rows.append({
                "species": s, "strain": j + 1, "strains": k, "depth": depth, "depth_lo": lo,
                "depth_hi": hi, "present_prob": present, "mixture_prob": mixture,
                **_placement(sp, one, st.points[:, j]),
                "survival_scale": float(st.diagnostics.get("survival_scale",
                                                           one.extra["survival_scale"])),
                "log_bf": one.log_bf, **checks, "check_failed": failed,
                **{key: v for key, v in st.diagnostics.items() if key != "survival_scale"},
            })  # fmt: skip
            grid = one.grid
            draws.append(pl.DataFrame({
                "species": np.full(len(lam), s, dtype=np.uint32),
                "strain": np.full(len(lam), j + 1, dtype=np.uint8),
                "draw": np.arange(len(lam), dtype=np.uint32),
                "ref_g": [sp.names[i] for i in grid.g[st.points[:, j]]],
                "ref_h": [sp.names[i] for i in grid.h[st.points[:, j]]],
                "t": grid.t[st.points[:, j]], "ell": grid.ell[st.points[:, j]], "lambda": lam,
            }))  # fmt: skip
            units.append(pl.DataFrame({
                "species": np.full(len(sp.units), s, dtype=np.uint32),
                "strain": np.full(len(sp.units), j + 1, dtype=np.uint8),
                "unit": sp.units.astype(np.uint32), "carriage_prob": st.carriage[j],
                "depth": st.depth[j], "survival": st.survival, "q": sp.q,
            }))  # fmt: skip
    table = pl.DataFrame(rows) if rows else pl.DataFrame(schema={"species": pl.UInt32})
    if rows:
        table = table.with_columns(
            pl.col("species").cast(pl.UInt32),
            relative_abundance=pl.col("depth") / pl.col("depth").sum(),
        )
        evidence = greedy_evidence(profile, Evidence.from_profile(profile), panel, detected)
        table = table.join(names, on="species").join(evidence, on="species", how="left")
        first = ["species", "id", "name", "strain", "strains", "depth", "depth_lo", "depth_hi",
                 "relative_abundance", "present_prob", "mixture_prob"]  # fmt: skip
        table = table.select(*first, pl.exclude(*first, "taxonomy"), "taxonomy").sort(
            "depth", descending=True
        )
    candidates = pl.DataFrame([
        {"species": s, "log_bf": f.log_bf, "present_prob": f.present_prob, "beta0": f.beta0,
         "beta1": f.beta1, **{k: v for k, v in f.extra.items() if k != "log_k0"}}
        for s, f in sorted(placed.fits.items())
    ]) if placed.fits else pl.DataFrame(schema={"species": pl.UInt32})  # fmt: skip
    if placed.fits:
        candidates = candidates.with_columns(pl.col("species").cast(pl.UInt32)).join(
            names, on="species"
        )
    return {
        "genome_profile": table,
        "candidates": candidates,
        "placements": pl.concat(draws) if draws else None,
        "strain_units": pl.concat(units) if units else None,
        "summary": placed.report | {"detected": len(detected), "strains": len(rows)},
    }
