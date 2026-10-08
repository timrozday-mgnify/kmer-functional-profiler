"""Evidence layer (phase 12, step 1): per-unit likelihoods rebuilt from the profile's columns,
as genotype likelihoods let imputation run downstream of calling.

A unit carried at depth *D* by a strain whose allele keeps a share *f* of the unit's *m* kept
k-mers: each kept k-mer is present with probability *f* (random survival, as the AAI model
has it) and, if present, hit ~ Poisson(*D*). With *k* hit k-mers and *H* hits on them,

    log L(D, f) = k log f + (m - k) log(1 - f (1 - e^-D)) + H log D - k D   (+ const)

(the binomial on hit k-mers times the zero-truncated Poisson on their hits; the log(1 - e^-D)
terms cancel). Its maximum is ``coverage_zi`` in *D* and ``present_zi`` in *f*. With *f*
profiled out it is the zero-truncated Poisson alone, whose likelihood-ratio interval is
``coverage_interval``. It is tempered by the profile's inflation (``coverage_zi_dispersion``
x read clumping), as ``coverage_interval`` is. *k* and *H* are the fit's: k = π m (1 - e^-c)
and H = π m c at ``coverage_zi`` c and ``present_zi`` π, so no index pass is needed.

- *Zero-hit units* (not in the profile): log L = m log(1 - f (1 - e^-D)) over the tier-2
  kept k-mers, exact under random survival (e^(-D m f) when f = 1).
- *Not carried:* a hit unit is then present by something outside the panel (prior β) at an
  unknown depth, with the generic marginal M = ∫ L dD df (log-uniform *D*, uniform *f*), or
  absent and hit by background. ``present_llr`` is log M / P(evidence | absent), so
  L(not carried) / M = β + (1 - β) e^-llr. A zero-hit unit not carried has no hits unless
  present from outside: 1 - β (1 - ε), ε = P(no hits | present) under the generic prior.

Everything a strain model needs is the ratio L(carried) / L(not carried) per unit
(:meth:`Evidence.log_ratio`): carriage is then Bernoulli mixing with 1.

Strain mixtures (opt-in sidecar, ``query --own-hist``): the hits of each unit's own k-mers,
a mixture of zero-truncated negative binomials at rates the strains fix, weights free
(:func:`mixture_loglik`); with one rate it is ``rate_mixture``'s one-rate fit.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl
from scipy.special import logsumexp
from scipy.stats import norm

from kmer_functional_profiler.query import (
    _ztnb_logpmf,
    mix_shape,
    ztp_lambda,
)

D_GRID: Final = np.geomspace(1e-2, 1e3, 81)  # the generic depth prior (log-uniform)
F_GRID: Final = (np.arange(40) + 0.5) / 40  # the generic survival prior (uniform)
F_FLOOR: Final = 1e-3  # survival below this is taken as this (an allele always keeps some)
CHUNK: Final = 2**21  # grid cells per batch of the generic marginal
# The background slab: units present from outside the panel at a depth log-normal around one
# of these (sd BG_SD in ln), the centre fitted per species by marginal likelihood (place.py).
BG_DEPTHS: Final = np.geomspace(0.01, 100, 17)
BG_SD: Final = 1.0


def log_likelihood(
    k: np.ndarray | float,
    h: np.ndarray | float,
    m: np.ndarray | float,
    d: np.ndarray | float,
    f: np.ndarray | float,
    w: np.ndarray | float = 1.0,
) -> np.ndarray:
    """log L(D, f) of a unit with ``k`` hit k-mers and ``h`` hits on ``m`` kept k-mers,
    tempered by ``w`` (module docstring); broadcasts. log(1 - f p) is computed as
    log((1 - f) + f e^-D), exact at f = 1."""
    d = np.maximum(np.asarray(d, dtype=np.float64), 1e-300)
    f = np.clip(np.asarray(f, dtype=np.float64), F_FLOOR, 1.0)
    k, h, m = (np.asarray(x, dtype=np.float64) for x in (k, h, m))
    with np.errstate(divide="ignore"):
        missed = np.logaddexp(np.log1p(-f), np.log(f) - d)
    term = np.where(k > 0, k * np.log(f) + h * np.log(d) - k * d, 0.0) + (m - k) * missed
    return np.asarray(term / np.asarray(w, dtype=np.float64))


def _background_log_prior() -> np.ndarray:
    """Per slab centre (``BG_DEPTHS``), the log prior weight of each ``D_GRID`` point:
    log-normal of sd ``BG_SD`` (in ln) x trapezoid on log D, normalised over the grid."""
    x = np.log(D_GRID)
    lp = norm.logpdf(x[None], np.log(BG_DEPTHS)[:, None], BG_SD) + np.log(np.gradient(x))
    return np.asarray(lp - logsumexp(lp, axis=1, keepdims=True))


def log_marginals(
    k: np.ndarray, h: np.ndarray, m: np.ndarray, w: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Per unit, log M (L averaged over the generic prior: log-uniform ``D_GRID`` x uniform
    ``F_GRID``) and log M_j under each background slab (log-normal *D* around
    ``BG_DEPTHS``[j], uniform *f*)."""
    gen, bg = np.empty(len(k)), np.empty((len(k), len(BG_DEPTHS)))
    prior = _background_log_prior()
    cells = len(D_GRID) * len(F_GRID)
    step = max(1, CHUNK // cells)
    d, f = D_GRID[None, :, None], F_GRID[None, None, :]
    for i in range(0, len(k), step):
        s = slice(i, i + step)
        ll = log_likelihood(k[s, None, None], h[s, None, None], m[s, None, None], d, f,
                            w[s, None, None])  # fmt: skip
        over_f = logsumexp(ll, axis=2) - np.log(len(F_GRID))  # units x D
        gen[s] = logsumexp(over_f, axis=1) - np.log(len(D_GRID))
        bg[s] = logsumexp(over_f[:, None, :] + prior[None], axis=2)
    return gen, bg


def log_marginal(k: np.ndarray, h: np.ndarray, m: np.ndarray, w: np.ndarray) -> np.ndarray:
    """log M per unit under the generic prior (:func:`log_marginals`)."""
    return log_marginals(k, h, m, w)[0]


def zero_hit_present(m: np.ndarray) -> np.ndarray:
    """ε: P(no hits | present) under the generic prior, for units of ``m`` kept k-mers."""
    return zero_hit_marginals(m)[0]


def zero_hit_marginals(m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """P(no hits | present) for units of ``m`` kept k-mers, under the generic prior and
    under each background slab (units x slabs)."""
    values, inverse = np.unique(np.asarray(m, dtype=np.float64), return_inverse=True)
    gen, bg = log_marginals(np.zeros(len(values)), np.zeros(len(values)), values,
                            np.ones(len(values)))  # fmt: skip
    return np.exp(gen)[inverse], np.exp(bg)[inverse]


def prior_odds(profile: pl.DataFrame) -> float:
    """The odds present_prob was fitted at: present_prob = expit(present_llr + log odds),
    one value per sample (median over units where both are finite and not at 0 or 1)."""
    rows = profile.filter(
        pl.col("present_llr").is_not_null(),
        pl.col("present_llr").abs() < 50,
        pl.col("present_prob").is_between(1e-12, 1 - 1e-12, closed="none"),
    )
    if rows.height == 0:
        return 1e-4
    p = rows["present_prob"].to_numpy()
    return float(np.exp(np.median(np.log(p / (1 - p)) - rows["present_llr"].to_numpy())))


@dataclass(frozen=True)
class Evidence:
    """The hit units of a profile as per-unit likelihood terms (module docstring).
    ``informative`` is False for units gather explained away (coverage 0): their own k-mers
    were not hit, but their k-mers are shared, so they get ratio 1 (no evidence) rather than
    a zero-hit term."""

    unit: np.ndarray  # index unit ids, sorted
    k: np.ndarray
    h: np.ndarray
    m: np.ndarray
    w: np.ndarray
    llr: np.ndarray
    present_prob: np.ndarray
    coverage: np.ndarray
    informative: np.ndarray
    odds: float  # the profile's prior odds of presence (:func:`prior_odds`)

    @classmethod
    def from_profile(cls, profile: pl.DataFrame) -> "Evidence":
        need = {"unit", "hits", "reads", "m_g", "coverage_zi", "present_zi", "present_prob",
                "present_llr"}  # fmt: skip
        missing = need - set(profile.columns)
        if missing:
            raise ValueError(f"the profile lacks columns {sorted(missing)} (query of phase 12)")
        prof = profile.filter(pl.col("hits") > 0).sort("unit")
        m_g = prof["m_g"].to_numpy().astype(np.float64)
        m = m_g
        if "m_dense" in prof.columns:  # the fit's tier
            dense = prof["m_dense"].fill_null(0).to_numpy().astype(np.float64)
            m = np.where(dense > 0, dense, m_g)
        c = prof["coverage_zi"].to_numpy().astype(np.float64)
        pi = prof["present_zi"].to_numpy().astype(np.float64)
        per_read = prof["hits"].to_numpy() / np.maximum(prof["reads"].to_numpy(), 1)
        clump = 1 + ztp_lambda(per_read.astype(np.float64)) * m / np.maximum(m_g, 1)
        dispersion = (
            prof["coverage_zi_dispersion"].fill_null(1.0).to_numpy()
            if "coverage_zi_dispersion" in prof.columns
            else np.ones(prof.height)
        )
        llr = prof["present_llr"].to_numpy().astype(np.float64)
        informative = (c > 0) & np.isfinite(llr)
        return cls(
            unit=prof["unit"].to_numpy().astype(np.int64),
            k=pi * m * -np.expm1(-c),
            h=pi * m * c,
            m=m,
            w=np.maximum(dispersion * clump, 1.0),
            llr=np.where(np.isfinite(llr), llr, -np.inf),
            present_prob=prof["present_prob"].to_numpy().astype(np.float64),
            coverage=c,
            informative=informative,
            odds=prior_odds(prof),
        )

    def index(self, units: np.ndarray) -> np.ndarray:
        """Row of each of ``units`` among the hit units, -1 where it had no hits."""
        units = np.asarray(units, dtype=np.int64)
        if len(self.unit) == 0:
            return np.full(len(units), -1)
        at = np.minimum(np.searchsorted(self.unit, units), len(self.unit) - 1)
        return np.where(self.unit[at] == units, at, -1)


@dataclass(frozen=True)
class UnitTerms:
    """The evidence of a fixed set of units (hit or not), ready for :meth:`log_ratio`:
    hit units' k, H, m, w, llr and log M (generic, and per background slab, relative to the
    generic); zero-hit units' tier-2 m and ε (generic, and per slab)."""

    hit: np.ndarray  # bool: unit is hit and informative
    blank: np.ndarray  # bool: hit but explained away (ratio 1)
    k: np.ndarray
    h: np.ndarray
    m: np.ndarray
    w: np.ndarray
    log_m: np.ndarray
    llr: np.ndarray
    eps: np.ndarray
    log_m_bg: np.ndarray  # units x slabs: log M_j - log M (hit units)
    eps_bg: np.ndarray  # units x slabs (zero-hit units)

    @classmethod
    def build(cls, ev: Evidence, units: np.ndarray, m_tier2: np.ndarray) -> "UnitTerms":
        """Terms of ``units`` (index ids) with ``m_tier2`` their tier-2 kept k-mers (from
        the panel's units table)."""
        at = ev.index(units)
        seen = at >= 0
        rows = np.where(seen, at, 0)
        hit = seen & ev.informative[rows]
        blank = seen & ~hit
        k = np.where(hit, ev.k[rows], 0.0)
        h = np.where(hit, ev.h[rows], 0.0)
        m = np.where(hit, ev.m[rows], np.asarray(m_tier2, dtype=np.float64))
        w = np.where(hit, ev.w[rows], 1.0)
        log_m = np.zeros(len(units))
        log_m_bg = np.zeros((len(units), len(BG_DEPTHS)))
        if hit.any():
            gen, bg = log_marginals(k[hit], h[hit], m[hit], w[hit])
            log_m[hit] = gen
            log_m_bg[hit] = bg - gen[:, None]
        zero = ~hit & ~blank
        eps = np.zeros(len(units))
        eps_bg = np.zeros((len(units), len(BG_DEPTHS)))
        if zero.any():
            eps[zero], eps_bg[zero] = zero_hit_marginals(m[zero])
        return cls(hit, blank, k, h, m, w, log_m, np.where(hit, ev.llr[rows], 0.0), eps,
                   log_m_bg, eps_bg)  # fmt: skip

    def log_not_carried(self, beta: float, slab: int | None = None) -> np.ndarray:
        """log L(not carried) per unit, hit units relative to M: present from outside the
        panel (prior ``beta``; depth from background ``slab``, else the generic prior), or
        absent."""
        rel = self.log_m_bg[:, slab] if slab is not None else 0.0
        eps = self.eps_bg[:, slab] if slab is not None else self.eps
        with np.errstate(divide="ignore"):
            hit = np.logaddexp(np.log(beta) + rel, np.log1p(-beta) - self.llr)
        return np.asarray(
            np.where(self.hit, hit, np.where(self.blank, 0.0, np.log1p(-beta * (1 - eps))))
        )

    def log_carried(
        self, d: np.ndarray, f: np.ndarray, rows: np.ndarray | slice = slice(None)
    ) -> np.ndarray:
        """log L(carried at depth ``d``, survival ``f``) of units ``rows`` (``d`` and ``f``
        broadcast against them on the first axis), hit units relative to M."""
        shape = (-1,) + (1,) * (np.ndim(d) - 1)
        k, h, m, w, log_m, blank = (
            x[rows].reshape(shape) for x in (self.k, self.h, self.m, self.w, self.log_m, self.blank)
        )
        return np.where(blank, 0.0, log_likelihood(k, h, m, d, f, w) - log_m)

    def log_ratio(
        self,
        d: np.ndarray,
        f: np.ndarray,
        beta: float,
        rows: np.ndarray | slice = slice(None),
        slab: int | None = None,
    ) -> np.ndarray:
        """log L(carried at ``d``, ``f``) / L(not carried) of units ``rows``."""
        shape = (-1,) + (1,) * (np.ndim(d) - 1)
        nc = self.log_not_carried(beta, slab)[rows].reshape(shape)
        return np.asarray(self.log_carried(d, f, rows) - nc)


# --- Strain mixtures: the own-k-mer histogram sidecar ---


@dataclass(frozen=True)
class Histogram:
    """Per unit, the hits of the k-mers only it holds (``query --own-hist``): rows ``row``
    (position in ``unit``), ``h`` hits and ``c`` k-mers with that many, and the sample's
    shape ``v`` (squared CV of uneven coverage, as ``rate_mixture`` reads it)."""

    unit: np.ndarray
    row: np.ndarray
    h: np.ndarray
    c: np.ndarray
    v: float

    @classmethod
    def read(cls, path: str | Path) -> "Histogram":
        return cls.from_frame(pl.read_parquet(path))

    @classmethod
    def from_frame(cls, table: pl.DataFrame) -> "Histogram":
        table = table.sort("unit", "hits")
        unit, row = np.unique(table["unit"].to_numpy().astype(np.int64), return_inverse=True)
        h = table["hits"].to_numpy().astype(np.float64)
        c = table["kmers"].to_numpy().astype(np.float64)
        return cls(unit, row.astype(np.int64), h, c, mix_shape(row, h, c) if len(row) else 0.0)


def mixture_loglik(
    hist: Histogram, rates: np.ndarray, iterations: int = 30, weights: np.ndarray | None = None
) -> np.ndarray:
    """Per unit of ``hist``, the log-likelihood of its own k-mers' hits as a mixture of
    zero-truncated negative binomials at fixed ``rates`` (units x r; equal rates are one
    component), the weights maximised by EM (or fixed, ``weights``). Without the log h!
    terms, as ``rate_mixture``'s."""
    n, r = rates.shape
    lp_c = _ztnb_logpmf(hist.h[:, None], np.maximum(rates[hist.row], 1e-9), hist.v)
    w = np.full((n, r), 1 / r) if weights is None else weights
    ll = np.zeros(n)
    for _ in range(1 if weights is not None else iterations):
        lp = np.log(np.maximum(w[hist.row], 1e-300)) + lp_c
        mx = lp.max(1, keepdims=True)
        resp = np.exp(lp - mx)
        s = resp.sum(1, keepdims=True)
        ll = np.bincount(hist.row, hist.c * (mx + np.log(s))[:, 0], n)
        if weights is None:
            resp /= s
            rs = np.stack([np.bincount(hist.row, hist.c * resp[:, j], n) for j in range(r)], 1)
            w = rs / np.maximum(rs.sum(1, keepdims=True), 1e-300)
    return ll
