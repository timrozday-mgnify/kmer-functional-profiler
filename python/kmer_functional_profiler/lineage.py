"""Continuous lineage model, prototype (phase 11, step 12; plan: Continuous lineage model).

A species' genomes differ in accessory content, and close strains share it. Each genome g
gets coordinates x_g in R^d, and carriage is a logistic factor model fitted on the species
index's per-genome carriage: logit P(g carries u) = b_u + w_u . x_g. A sample holds K
lineages of the species, each a point x_k near the reference genomes (a Gaussian mixture
centred on them) at depth λ_k, carrying u with probability σ(b_u + w_u . x_k). The hits are
the species model's: h_u ~ Poisson(others_u + Σ_k z_{k,u} λ_k e_u), with the background free
to add hits (ℓ, :func:`~kmer_functional_profiler.species._explain`) and every likelihood
tempered by the species fit's φ. The other species are fixed at the hits the species fit
gave them (block coordinate ascent). Carriage is summed out over its 2^K patterns per unit;
λ and x are fitted by MAP (NumPyro SVI, ``AutoDelta``) from several starts.

It reduces to the species model at d = 0, K = 1. Needs the optional ``phylo`` dependency
group (JAX, NumPyro); nothing else in the package imports it.
"""

from dataclasses import dataclass
from itertools import product
from typing import Any, Final

import jax
import jax.numpy as jnp
import numpy as np
import numpyro
import numpyro.distributions as dist
import polars as pl
from jax.nn import log_sigmoid
from jax.scipy.special import gammaln, logsumexp
from numpyro.infer import SVI, Trace_ELBO
from numpyro.infer.autoguide import AutoDelta
from numpyro.infer.initialization import init_to_value

from kmer_functional_profiler.species import (
    BG_PRIOR,
    CORE,
    DETECTED,
    PRIOR_PRESENT,
    SpeciesIndex,
    _log_nb,
    fit_species,
)

jax.config.update("jax_enable_x64", True)  # type: ignore[no-untyped-call]

# ponytail: guesses for a prototype, to tune on the strain hold-out benchmark (step 10)
DIM: Final = 2  # lineage coordinates
LINEAGES: Final = 2  # K per species
STARTS: Final = 4  # MAP starts per species, each at K random reference genomes
STEPS: Final = 1500  # Adam steps per fit
LEARNING_RATE: Final = 0.05
DEPTH_SHAPE: Final = 0.5  # Gamma prior on a lineage's depth: favours few lineages
DEPTH_RATE: Final = 1e-3
MIN_SHARE: Final = 0.1  # lineages below this share of their species' depth are dropped
NEAREST: Final = 3  # reference genomes reported per lineage
LINEAGE_COLUMNS: Final = ("species", "id", "name", "lineage", "depth", "share", "nearest",
                           "distance", "accessory_called")  # fmt: skip
UNIT_COLUMNS: Final = ("species", "lineage", "unit", "prevalence", "hits", "carriage_prob",
                        "hits_assigned")  # fmt: skip


@dataclass
class Factors:
    """A species' logistic factor model: ``b`` (U), ``w`` (U x d), genome coordinates
    ``x`` (G x d) and ``sigma``, the spread of a lineage around a reference genome (the
    median distance from a genome to its nearest neighbour)."""

    b: np.ndarray
    w: np.ndarray
    x: np.ndarray
    sigma: float


def _factor_model(carried: jnp.ndarray, dim: int) -> None:
    g, u = carried.shape
    b = numpyro.sample("b", dist.Normal(0.0, 3.0).expand([u]).to_event(1))
    w = numpyro.sample("w", dist.Normal(0.0, 1.0).expand([u, dim]).to_event(2))
    x = numpyro.sample("x", dist.Normal(0.0, 1.0).expand([g, dim]).to_event(2))
    numpyro.sample("c", dist.Bernoulli(logits=b + x @ w.T).to_event(2), obs=carried)


def _map(model: Any, init: dict[str, Any], seed: int, steps: int, *args: Any) -> tuple[Any, float]:
    """MAP point and final loss of ``model`` by SVI with an ``AutoDelta`` guide."""
    guide = AutoDelta(model, init_loc_fn=init_to_value(values=init))
    svi = SVI(model, guide, numpyro.optim.Adam(LEARNING_RATE), Trace_ELBO())
    result = svi.run(jax.random.PRNGKey(seed), steps, *args, progress_bar=False)
    point = guide.median(result.params)
    return point, float(result.losses[-1])


def fit_factors(carried: np.ndarray, dim: int = DIM, seed: int = 0, steps: int = STEPS) -> Factors:
    """The logistic factor model of a (genomes x units) 0/1 carriage matrix, by MAP. ``b``
    starts at each unit's logit frequency, ``w`` and ``x`` at random (zero is a saddle)."""
    freq = np.clip(carried.mean(axis=0), 0.01, 0.99)
    init = {"b": jnp.asarray(np.log(freq / (1 - freq)))}
    point, _ = _map(_factor_model, init, seed, steps, jnp.asarray(carried, float), dim)
    x = np.asarray(point["x"])
    distance = np.linalg.norm(x[:, None] - x[None], axis=-1)
    np.fill_diagonal(distance, np.inf)
    sigma = float(np.median(distance.min(axis=1))) if len(x) > 1 else 1.0
    return Factors(np.asarray(point["b"]), np.asarray(point["w"]), x, max(sigma, 0.1))


def _patterns(k: int) -> jnp.ndarray:
    """Every carriage pattern of k lineages, (2^k x k)."""
    return jnp.asarray(list(product((0.0, 1.0), repeat=k)))


def _log_pois(h: jnp.ndarray, m: jnp.ndarray) -> jnp.ndarray:
    return h * jnp.log(jnp.maximum(m, 1e-12)) - m - gammaln(h + 1)


def _joint(
    lam: jnp.ndarray, x: jnp.ndarray, data: dict[str, Any], f: Factors
) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Per (pattern, unit): log prior of the pattern plus ℓ of the unit's hits given it;
    and the pattern matrix and the lineages' means, (2^K x U) and (K x U)."""
    z = _patterns(len(lam))
    logits = jnp.asarray(f.b) + x @ jnp.asarray(f.w).T  # K x U
    prior = z @ log_sigmoid(logits) + (1 - z) @ log_sigmoid(-logits)
    own = lam[:, None] * data["e"][None]
    mean = data["others"][None] + z @ own
    h, phi, log_bg = data["h"], data["phi"], data["log_bg"]
    ell = jnp.logaddexp(
        jnp.log1p(-jnp.exp(log_bg)) + _log_pois(h, mean) / phi,
        log_bg + (data["nb"] - mean) / phi,
    )
    return prior + ell, z, own


def _lineage_model(data: dict[str, Any], f: Factors, k: int) -> None:
    lam = numpyro.sample("lam", dist.Gamma(DEPTH_SHAPE, DEPTH_RATE).expand([k]).to_event(1))
    x = numpyro.sample("x", dist.ImproperUniform(dist.constraints.real, (), (k, f.x.shape[1])))
    # each lineage near some reference genome: a mixture of N(x_g, σ²) over the genomes
    near = dist.Normal(jnp.asarray(f.x), f.sigma).log_prob(x[:, None]).sum(-1)  # K x G
    numpyro.factor("near", (logsumexp(near, axis=1) - jnp.log(f.x.shape[0])).sum())
    joint, _, _ = _joint(lam, x, data, f)
    numpyro.factor("hits", logsumexp(joint, axis=0).sum())


def fit_lineages(
    h: np.ndarray,
    others: np.ndarray,
    e: np.ndarray,
    f: Factors,
    *,
    phi: float = 1.0,
    background_prior: float = BG_PRIOR,
    k: int = LINEAGES,
    starts: int = STARTS,
    steps: int = STEPS,
    seed: int = 0,
) -> dict[str, np.ndarray]:
    """K lineages of a species from its units' hits ``h``, the other sources' hits
    ``others`` and its expected hits per genome copy ``e`` (all over the factor model's
    units). Returns ``depth`` (K), ``x`` (K x d), ``carriage`` (K x U: P(z = 1 | h)) and
    ``assigned`` (K x U: expected hits), lineages by decreasing depth."""
    with np.errstate(divide="ignore"):
        log_bg = float(np.log(background_prior))
    data: dict[str, Any] = {"h": jnp.asarray(h), "others": jnp.asarray(others), "e": jnp.asarray(e),
            "nb": jnp.asarray(_log_nb(h)), "phi": phi, "log_bg": log_bg}  # fmt: skip
    rng = np.random.default_rng(seed)
    exposure = float((e * (1 / (1 + np.exp(-f.b)))).sum())
    depth0 = max(float(np.maximum(h - others, 0).sum()) / max(exposure, 1e-12), 1e-3)
    best: tuple[Any, float] | None = None
    for start in range(starts):
        picks = rng.choice(len(f.x), size=k, replace=len(f.x) < k)
        init = {"lam": jnp.full(k, depth0 / k), "x": jnp.asarray(f.x[picks])}
        point, loss = _map(_lineage_model, init, seed + start, steps, data, f, k)
        if best is None or loss < best[1]:
            best = (point, loss)
    assert best is not None
    lam, x = best[0]["lam"], best[0]["x"]
    joint, z, own = _joint(lam, x, data, f)
    post = jnp.exp(joint - logsumexp(joint, axis=0))  # 2^K x U
    carriage = z.T @ post
    total = jnp.maximum(data["others"] + z @ own, 1e-12)
    share = z[:, :, None] * own[None] / total[:, None]  # 2^K x K x U
    assigned = (post[:, None] * share).sum(0) * data["h"][None]
    order = np.argsort(-np.asarray(lam))
    fit = {"depth": lam, "x": x, "carriage": carriage, "assigned": assigned}
    return {key: np.asarray(value)[order] for key, value in fit.items()}


def _prune(fit: dict[str, np.ndarray], sigma: float) -> dict[str, np.ndarray]:
    """Lineages within σ of a deeper one merge into it (depth and hits summed, its carriage
    kept); then lineages under ``MIN_SHARE`` of the species' depth are dropped."""
    keep: list[int] = []
    depth, assigned = fit["depth"].copy(), fit["assigned"].copy()
    for i in range(len(depth)):
        near = [j for j in keep if np.linalg.norm(fit["x"][i] - fit["x"][j]) < sigma]
        if near:
            depth[near[0]] += depth[i]
            assigned[near[0]] += assigned[i]
        else:
            keep.append(i)
    keep = [i for i in keep if depth[i] >= MIN_SHARE * depth.sum()]
    return {"depth": depth[keep], "x": fit["x"][keep], "carriage": fit["carriage"][keep],
            "assigned": assigned[keep]}  # fmt: skip


def _carriage(si: SpeciesIndex, species_id: str, units: np.ndarray) -> tuple[list[str], np.ndarray]:
    """Names of a species' genomes and their (genomes x ``units``) carriage."""
    genomes = pl.read_parquet(si.path / "genomes.parquet").filter(pl.col("species") == species_id)
    rows = (
        pl.scan_parquet(si.path / "genome_units.parquet")
        .filter(pl.col("genome").is_in(genomes["genome"].implode()))
        .collect()
    )
    g = np.searchsorted(genomes["genome"].to_numpy(), rows["genome"].to_numpy())
    unit = rows["unit"].to_numpy()
    column = np.searchsorted(units, unit)
    ok = (column < len(units)) & (units[np.minimum(column, len(units) - 1)] == unit)
    carried = np.zeros((genomes.height, len(units)))
    carried[g[ok], column[ok]] = 1.0
    return genomes["name"].to_list(), carried


def lineage_profile(
    prof: pl.DataFrame,
    si: SpeciesIndex,
    *,
    dim: int = DIM,
    k: int = LINEAGES,
    starts: int = STARTS,
    steps: int = STEPS,
    seed: int = 0,
    prior: float = PRIOR_PRESENT,
    background_prior: float = BG_PRIOR,
    **screen: Any,
) -> dict[str, Any]:
    """The species fit, then lineages within each detected species with ≥ 2 genomes
    (one-genome species keep one lineage at the species fit's depth and carriage).

    Returns ``lineages`` (``species``, ``id``, ``name``, ``lineage``: 0 the deepest,
    ``depth``, ``share`` of the species' depth, ``nearest``: the ``NEAREST`` closest
    reference genomes, ``distance`` to the closest in σ, ``accessory_called``), ``units``
    (``species``, ``lineage``, ``unit``, ``prevalence``, ``hits``, ``carriage_prob``,
    ``hits_assigned``) and ``summary``."""
    species, pairs, report = fit_species(
        prof, si, prior=prior, background_prior=background_prior, **screen
    )
    phi = float(report["fit_phi"])
    detected = species.filter(pl.col("present_prob") >= DETECTED)
    names = si.species.select(pl.col("species").cast(pl.UInt32), "id", "name")
    lineages, units = [], []
    for s, depth in detected.select("species", "depth").iter_rows():
        sid, name = names.filter(pl.col("species") == s).row(0)[1:]
        own = pairs.filter(pl.col("species") == s).sort("unit")
        unit = own["unit"].to_numpy()
        genomes, carried = _carriage(si, sid, unit)
        h, q = own["hits"].to_numpy(), own["q"].to_numpy()
        if len(genomes) < 2:
            fit = {"depth": np.array([depth]), "carriage": own["carriage_prob"].to_numpy()[None],
                   "assigned": own["hits_assigned"].to_numpy()[None]}  # fmt: skip
            nearest, distance = [",".join(genomes)], [0.0]
        else:
            f = fit_factors(carried, dim, seed, steps)
            others = np.maximum(h - own["hits_assigned"].to_numpy(), 0.0)
            fitted = fit_lineages(
                h, others, own["e"].to_numpy(), f, phi=phi, background_prior=background_prior,
                k=k, starts=starts, steps=steps, seed=seed,
            )  # fmt: skip
            fit = _prune(fitted, f.sigma)
            d = np.linalg.norm(fit["x"][:, None] - f.x[None], axis=-1) / f.sigma
            nearest = [",".join(genomes[j] for j in np.argsort(row)[:NEAREST]) for row in d]
            distance = list(d.min(axis=1))
        n = len(fit["depth"])
        lineages.append(pl.DataFrame({
            "species": [s] * n, "id": [sid] * n, "name": [name] * n,
            "lineage": np.arange(n, dtype=np.uint32), "depth": fit["depth"],
            "share": fit["depth"] / fit["depth"].sum(), "nearest": nearest, "distance": distance,
            "accessory_called": ((fit["carriage"] >= 0.5) & (q < CORE)).sum(axis=1),
        }))  # fmt: skip
        units.append(pl.DataFrame({
            "species": np.repeat(s, n * len(unit)).astype(np.uint32),
            "lineage": np.repeat(np.arange(n, dtype=np.uint32), len(unit)),
            "unit": np.tile(unit, n), "prevalence": np.tile(q, n), "hits": np.tile(h, n),
            "carriage_prob": fit["carriage"].ravel(), "hits_assigned": fit["assigned"].ravel(),
        }))  # fmt: skip
    if not lineages:  # nothing detected: headers only
        lineages = [pl.DataFrame(schema=dict.fromkeys(LINEAGE_COLUMNS, pl.Float64))]
        units = [pl.DataFrame(schema=dict.fromkeys(UNIT_COLUMNS, pl.Float64))]
    table = pl.concat(lineages).with_columns(
        relative_abundance=pl.col("depth") / pl.col("depth").sum()
    )
    return {
        "lineages": table,
        "units": pl.concat(units),
        "summary": report | {"species_detected": detected.height, "lineages": table.height},
    }
