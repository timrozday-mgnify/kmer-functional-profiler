"""Strain mixtures (phase 12, step 4): K ≥ 2 strains of a species by MCMC, K by
stepping-stone sampling.

The state is each strain's placement θ_k (a point of the K = 1 grid) and depth λ_k. Carriage
is summed out exactly per unit: each strain is in one of its states (not carrying the unit,
or carrying it copied from g, h or the species average, weights as for one strain), and the
unit's likelihood sums over every combination. Where several strains carry a unit its k-mers
are in any of their alleles, at the depth that keeps the expected hits (as for species
sharing units). With the own-k-mer histogram (``query --own-hist``), a unit carried by
several strains also gets its histogram's likelihood as a mixture at the rates of every
non-empty subset of them, against one rate at their summed depth.

Moves: slice sampling on each log λ_k; Metropolis-Hastings on θ_k, either within its
neighbour pair (another t or ℓ, uniformly) or by an independence proposal from the K = 1
posterior, the K = 1 posterior of the residual (a strain fitted with the first at its
posterior mean) and the prior. Strains are exchangeable; samples are reported sorted by
depth.

The marginal likelihood of K strains is by generalised stepping-stone sampling (Fan et al.
2011): the path p_β ∝ (L prior)^β r^(1 - β) runs from a reference r close to the posterior
(the K = 1 and residual fits' placements, log λ normal around their depths, symmetrised over
labels), sampled exactly, to the posterior. (From the prior, as first written, the first
rung's draws were hopeless and K ≥ 2 came out ~1000 nats low.) For K = 1 it is checked
against the exact grid. R-hat and effective sample sizes of the log-likelihood and each
sorted log λ are reported.
"""

import itertools
from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Final

import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm

from kfp_genomes.evidence import Histogram, UnitTerms
from kfp_genomes.panel import SpeciesPanel
from kfp_genomes.place import (
    COPY_CHANGE,
    LAMBDA_RANGE,
    LAMBDA_SD,
    SCALE_GRID,
    OneStrain,
    Others,
    Placements,
    _not_carried,
    copy_change,
    fit_one,
    survival,
)
from kmer_functional_profiler.query import MIX_SHAPES, _ztnb_logpmf, ztp_lambda

TEMPERATURES: Final = 8  # stepping-stone rungs (β = (j / J)^3)
SWEEPS: Final = 30  # sweeps per rung (after BURN)
BURN: Final = 5
CHAINS: Final = 4
SLICE_WIDTH: Final = 0.03  # on log λ at β = 1 (plus the reference's sd x (1 - β))
SLICE_STEPS: Final = 20
PROPOSAL_PRIOR: Final = 0.1  # independence proposal: share of the prior
PLACEMENT_MOVES: Final = 3  # per strain and sweep: one within its pair, the rest independent
SEARCH_ROUNDS: Final = 4  # coordinate-ascent rounds for the K ≥ 2 starts
SPLITS: Final = (0.35, 0.5, 0.65, 0.8, 0.9)  # shares of a strain's depth tried as a split
SCALE_DROPS: Final = np.array([0.05, 0.1, 0.15, 0.2])  # survival scales tried with a split
SEARCH_DEPTHS: Final = 9  # λ points per search (a coarse log grid, then a finer one)
PAIR_TOPS: Final = 3  # reference pairs of the K = 1 posterior tried as two strains
PAIR_ELL: Final = 0.25
SEARCH_ELL: Final = (0.0, 0.25, 0.5)  # pulls of the nodes the search tries (new strains
# differ from every reference: at ℓ = 0 a node forbids the units it lacks)
PILOT: Final = 20  # sweeps per chain of the K ≥ 2 survival-scale pilot
REFERENCE_SD: Final = 0.3  # the stepping-stone reference's sd of log λ around each start
SAMPLES: Final = 200  # posterior samples kept (pooled, thinned) for summaries
STRAIN_RATIO: Final = 0.2  # P(K = k + 1 | K > k): truncated geometric prior on K ≥ 1


@dataclass(frozen=True)
class _States:
    """One strain's states on every unit: log weight of not carrying, and per source the
    log weight of carrying from it, with its copies and survival."""

    log_not: np.ndarray
    log_w: np.ndarray  # sources x U
    n: np.ndarray
    f: np.ndarray


class Model:
    """The K-strain likelihood of one species (module docstring)."""

    def __init__(
        self,
        sp: SpeciesPanel,
        terms: UnitTerms,
        grid: Placements,
        not_carried: np.ndarray,
        hist: Histogram | None = None,
        scale: float = 1.0,
    ) -> None:
        self.sp, self.terms, self.grid, self.nc = sp, terms, grid, not_carried
        self.rows = np.arange(len(sp.units))
        self.xt = np.vstack([sp.xt, sp.q[None]])
        self.n_src = np.vstack([sp.n, sp.n_mean[None]])
        self.scale = scale
        self.f_src = survival(sp, scale)
        self._cache: dict[int, _States] = {}
        self.hist = None
        if hist is not None:  # the histogram rows of this species' units
            at = np.searchsorted(sp.units, hist.unit)
            ok = (at < len(sp.units)) & (sp.units[np.minimum(at, len(sp.units) - 1)] == hist.unit)
            keep = ok[hist.row]
            if keep.any():
                # rows renumbered to the kept units, which map to this species' units
                kept_rows = np.flatnonzero(ok)
                row = np.searchsorted(kept_rows, hist.row[keep])
                self.hist = (at[ok], row, hist.h[keep], hist.c[keep], hist.v)

    def set_shape(self, v: float) -> None:
        """Another histogram shape (squared CV of one allele's uneven coverage)."""
        if self.hist is not None:
            self.hist = (*self.hist[:4], v)

    def set_scale(self, scale: float) -> None:
        """Another survival scale (the strains' alleles' survival relative to the panel's)."""
        self.scale = scale
        self.f_src = survival(self.sp, scale)
        self._cache.clear()

    def states(self, point: int) -> _States:
        if point not in self._cache:
            wg, wh, wa = (float(x[0]) for x in self.grid.weights(np.array([point])))
            g, h, avg = int(self.grid.g[point]), int(self.grid.h[point]), len(self.n_src) - 1
            src = [(g, wg + (wh if h == g else 0.0))] + ([(h, wh)] if h != g else [])
            src = [(s, w) for s, w in src + [(avg, wa)] if w > 0]
            not_c = sum(w * (1 - self.xt[s]) for s, w in src)
            with np.errstate(divide="ignore"):
                self._cache[point] = _States(
                    np.log(not_c),
                    np.log(np.stack([w * self.xt[s] for s, w in src])),
                    np.stack([self.n_src[s] for s, _ in src]),
                    np.stack([self.f_src[s] for s, _ in src]),
                )
        return self._cache[point]

    def _combos(self, states: list[_States]) -> Iterator[tuple[int, ...]]:
        """Every combination of per-strain states (-1: not carried, else a source)."""
        return itertools.product(*[range(-1, len(s.log_w)) for s in states])

    def unit_terms(
        self, points: list[int], lam: np.ndarray
    ) -> tuple[np.ndarray, list[tuple[int, ...]]]:
        """log (weight x likelihood) per combination and unit, and the combinations."""
        states = [self.states(p) for p in points]
        combos = list(self._combos(states))
        out = np.empty((len(combos), len(self.rows)))
        f_mean = self.f_src[-1]
        self._change = [copy_change(self.terms, np.array([x]), f_mean)[:, 0] for x in lam]
        changes: dict[tuple[int, ...], np.ndarray] = {}  # per set of carrying strains
        for c, combo in enumerate(combos):
            log_w = np.zeros(len(self.rows))
            carriers: list[tuple[np.ndarray, np.ndarray, int]] = []
            for k, state in enumerate(combo):
                s = states[k]
                if state < 0:
                    log_w += s.log_not
                else:
                    log_w += s.log_w[state]
                    carriers.append((lam[k] * s.n[state], s.f[state], k))
            if not carriers:
                out[c] = log_w + self.nc
                continue
            if len(carriers) == 1:  # as one strain's, copy changes included
                ll = np.logaddexp(
                    np.log1p(-COPY_CHANGE)
                    + self.terms.log_carried(carriers[0][0], carriers[0][1], self.rows),
                    self._change[carriers[0][2]],
                )
            else:
                ll = self.terms.log_carried_alleles(
                    [c[0] for c in carriers], [c[1] for c in carriers], self.rows
                )
                if self.hist is not None:
                    ll = ll + self._histogram(carriers)
                # copy changes as for one strain, at the carriers' summed depth (without
                # them a strain at no depth "carrying" every unit escaped their cost)
                who = tuple(c[2] for c in carriers)
                if who not in changes:
                    total = float(sum(lam[k] for k in who))
                    changes[who] = copy_change(self.terms, np.array([total]), self.f_src[-1])[:, 0]
                ll = np.logaddexp(np.log1p(-COPY_CHANGE) + ll, changes[who])
            out[c] = log_w + ll
        return out, combos

    def _histogram(self, carriers: list[tuple[np.ndarray, np.ndarray, int]]) -> np.ndarray:
        """Per unit, the own-k-mer histogram's log-likelihood at the rates of every
        non-empty subset of the carriers (weights: the share of k-mers in exactly those
        alleles, from their survival, that are hit), less its one-rate likelihood at the
        rate the unit's likelihood takes
        for them (:func:`~kfp_genomes.evidence.log_likelihood_alleles`'s D*)."""
        assert self.hist is not None
        units, row, h, c, v = self.hist
        n = len(units)
        subsets = [s for r in range(1, len(carriers) + 1)
                   for s in itertools.combinations(range(len(carriers)), r)]  # fmt: skip
        rates = np.stack([sum(carriers[k][0][units] for k in s) for s in subsets], 1)
        f = np.stack([cf[1][units] for cf in carriers], 1)
        weights = np.stack(
            [np.prod([f[:, k] if k in s else 1 - f[:, k] for k in range(len(carriers))], 0)
             for s in subsets], 1,
        )  # fmt: skip
        # the histogram holds hit k-mers only: each subset's share given a hit
        p0 = np.exp(-rates) if v == 0 else (1 + v * rates) ** (-1 / v)
        weights = weights * (1 - p0)
        weights /= np.maximum(weights.sum(1, keepdims=True), 1e-300)
        lp = np.log(np.maximum(weights[row], 1e-300)) + _ztnb_logpmf(
            h[:, None], np.maximum(rates[row], 1e-9), v
        )
        p_hit = 1 - np.prod(1 - f * -np.expm1(-np.stack([cd[0][units] for cd in carriers], 1)), 1)
        mean = sum(cd[0][units] * cd[1][units] for cd in carriers)
        d_star = ztp_lambda(mean / np.maximum(p_hit, 1e-300))
        mix = np.bincount(row, c * logsumexp(lp, axis=1), n)
        one = np.bincount(row, c * _ztnb_logpmf(h, np.maximum(d_star[row], 1e-9), v), n)
        out = np.zeros(len(self.rows))
        np.add.at(out, units, mix - one)
        return out

    def loglik(self, points: list[int], lam: np.ndarray) -> float:
        terms, _ = self.unit_terms(points, lam)
        return float(logsumexp(terms, axis=0).sum())

    def carriage(self, points: list[int], lam: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """P(strain k carries unit u | θ, λ, data), and E[λ_k n | carries, θ, λ, data]
        times that probability: each strains x units."""
        terms, combos = self.unit_terms(points, lam)
        post = np.exp(terms - logsumexp(terms, axis=0, keepdims=True))
        states = [self.states(p) for p in points]
        p = np.zeros((len(points), len(self.rows)))
        depth = np.zeros_like(p)
        for c, combo in enumerate(combos):
            for k, state in enumerate(combo):
                if state >= 0:
                    p[k] += post[c]
                    depth[k] += post[c] * lam[k] * states[k].n[state]
        return p, depth


_LO, _HI = np.log(LAMBDA_RANGE)
_LOG_MASS = float(np.log(norm.cdf(_HI / LAMBDA_SD) - norm.cdf(_LO / LAMBDA_SD)))
_LOG_NORM = -0.5 * np.log(2 * np.pi) - np.log(LAMBDA_SD) - _LOG_MASS


def _log_lambda_prior(x: np.ndarray) -> np.ndarray:
    """log density of log λ under its prior (normal, truncated to ``LAMBDA_RANGE``)."""
    x = np.asarray(x, dtype=np.float64)
    inside = (x >= _LO) & (x <= _HI)
    return np.where(inside, _LOG_NORM - 0.5 * (x / LAMBDA_SD) ** 2, -np.inf)


@dataclass
class Chain:
    points: list[int]
    log_lam: np.ndarray
    ll: float


@dataclass
class MixtureFit:
    """K strains of a species: the log Bayes factor against K = 0, posterior samples sorted
    by depth, per-strain carriage, and diagnostics."""

    k: int
    log_bf: float
    points: np.ndarray  # samples x K
    lam: np.ndarray  # samples x K
    carriage: np.ndarray  # K x U: P(strain k carries u | data), Rao-Blackwellised
    depth: np.ndarray  # K x U: E[λ_k n | strain k carries u, data]
    diagnostics: dict[str, float] = field(default_factory=dict)


class Reference:
    """The stepping-stone reference: each strain's placement from the proposal ``q`` and
    log λ normal around a start, symmetrised over the strains' labels (the posterior is
    exchangeable). Sampled exactly; its density is the path's other end."""

    def __init__(self, q: np.ndarray, centres: np.ndarray, sd: float = REFERENCE_SD) -> None:
        self.q, self.log_q = q, np.log(q)
        self.centres, self.sd = centres, sd
        k = len(centres)
        self.perms = np.array(list(itertools.permutations(range(k))))
        self.log_k_factorial = float(np.log(len(self.perms)))

    def log_density(self, points: list[int], log_lam: np.ndarray) -> float:
        pts = np.asarray(points)[self.perms]  # perms x K
        lam = log_lam[self.perms]
        z = (lam - self.centres[None]) / self.sd
        terms = self.log_q[pts] - 0.5 * z**2 - np.log(self.sd) - 0.5 * np.log(2 * np.pi)
        return float(logsumexp(terms.sum(1)) - self.log_k_factorial)

    def draw(self, rng: np.random.Generator) -> tuple[list[int], np.ndarray]:
        k = len(self.centres)
        points = rng.choice(len(self.q), k, p=self.q)
        log_lam = rng.normal(self.centres, self.sd)
        order = rng.permutation(k)
        return [int(x) for x in points[order]], log_lam[order]


class Sampler:
    """MCMC for K strains of one species on the path p_β ∝ (L prior)^β r^(1 - β) between
    the reference r (β = 0) and the posterior (β = 1); module docstring."""

    def __init__(self, model: Model, k: int, ref: Reference, rng: np.random.Generator):
        self.model, self.k, self.ref, self.rng = model, k, ref, rng
        grid = model.grid
        key = grid.g * (grid.g.max() + 1) + grid.h  # the points of one neighbour pair
        self.groups = {int(u): np.flatnonzero(key == u) for u in np.unique(key)}
        self.group_of = key

    def log_prior(self, points: list[int], log_lam: np.ndarray) -> float:
        return float(self.model.grid.log_prior[points].sum() + _log_lambda_prior(log_lam).sum())

    def log_target(self, points: list[int], log_lam: np.ndarray, ll: float, beta: float) -> float:
        prior = self.log_prior(points, log_lam)
        if not np.isfinite(prior):
            return -np.inf
        ref = self.ref.log_density(points, log_lam) if beta < 1 else 0.0
        return beta * (ll + prior) + (1 - beta) * ref

    def _try(self, c: Chain, points: list[int], log_lam: np.ndarray, beta: float,
             log_ratio: float = 0.0) -> None:  # fmt: skip
        """Metropolis-Hastings to (points, log_lam), ``log_ratio`` the proposal's."""
        if not np.isfinite(_log_lambda_prior(log_lam).sum()):
            return
        ll = self.model.loglik(points, np.exp(log_lam))
        new = self.log_target(points, log_lam, ll, beta)
        old = self.log_target(c.points, c.log_lam, c.ll, beta)
        if np.log(self.rng.random()) < new - old + log_ratio:
            c.points, c.log_lam, c.ll = points, log_lam, ll

    def sweep(self, c: Chain, beta: float) -> Chain:
        for j in range(self.k):
            self._slice(c, j, beta)
            for move in range(PLACEMENT_MOVES):
                points = list(c.points)
                p_old = points[j]
                if move == 0:  # within the neighbour pair: another t or ℓ
                    points[j] = int(self.rng.choice(self.groups[int(self.group_of[p_old])]))
                    log_ratio = 0.0
                else:  # independence proposal
                    points[j] = int(self.rng.choice(len(self.ref.q), p=self.ref.q))
                    log_ratio = float(self.ref.log_q[p_old] - self.ref.log_q[points[j]])
                if points[j] != p_old:
                    self._try(c, points, c.log_lam.copy(), beta, log_ratio)
        if self.k > 1:  # swap two strains' depths
            a, b = self.rng.choice(self.k, 2, replace=False)
            log_lam = c.log_lam.copy()
            log_lam[[a, b]] = log_lam[[b, a]]
            self._try(c, list(c.points), log_lam, beta)
        return c

    def _slice(self, c: Chain, j: int, beta: float) -> None:
        """Slice sampling on log λ_j (stepping out)."""

        def logp(x: float) -> tuple[float, float]:
            log_lam = c.log_lam.copy()
            log_lam[j] = x
            if not np.isfinite(_log_lambda_prior(log_lam).sum()):
                return -np.inf, -np.inf
            ll = self.model.loglik(c.points, np.exp(log_lam))
            return self.log_target(c.points, log_lam, ll, beta), ll

        x0 = c.log_lam[j]
        width = REFERENCE_SD * (1 - beta) + SLICE_WIDTH  # the posterior is narrow at β = 1
        y = self.log_target(c.points, c.log_lam, c.ll, beta) - self.rng.exponential()
        left = x0 - width * self.rng.random()
        right = left + width
        for _ in range(SLICE_STEPS):
            if logp(left)[0] <= y:
                break
            left -= width
        for _ in range(SLICE_STEPS):
            if logp(right)[0] <= y:
                break
            right += width
        for _ in range(100):
            x1 = self.rng.uniform(left, right)
            value, ll = logp(x1)
            if value > y:
                c.log_lam[j], c.ll = x1, ll
                return
            if x1 < x0:
                left = x1
            else:
                right = x1


def _best_depth(model: Model, points: list[int], lam: np.ndarray, j: int) -> tuple[float, float]:
    """λ_j maximising the likelihood, the others fixed: a log grid over the range, then a
    finer one around its best; returns (λ_j, log-likelihood)."""
    lo, hi = np.log(LAMBDA_RANGE)
    best = (-np.inf, 0.0)
    for grid in (np.linspace(lo, hi, SEARCH_DEPTHS), None):
        if grid is None:
            step = (hi - lo) / (SEARCH_DEPTHS - 1)
            grid = np.clip(best[1] + np.linspace(-step, step, SEARCH_DEPTHS), lo, hi)
        for x in grid:
            trial = lam.copy()
            trial[j] = np.exp(x)
            ll = model.loglik(points, trial)
            if ll > best[0]:
                best = (ll, float(x))
    return float(np.exp(best[1])), best[0]


def _pair_splits(
    model: Model, one: OneStrain, pts: list[int], depth: np.ndarray, fit_scale: bool
) -> tuple[list[int], np.ndarray, float]:
    """One strain placed between two references (the K = 1 posterior's top points with
    g ≠ h) split into two, one at each reference (ℓ = ``PAIR_ELL``), with the survival
    scale lowered and each depth maximised in turn (the K = 1 depth is no guide: its
    raised survival absorbed part of the second strain); the best start found (points,
    depths, log-likelihood), else the given one."""
    grid = model.grid
    best = (model.loglik(pts, depth), list(pts), depth.copy(), model.scale)
    post = one.point_posterior()
    tops = [int(p) for p in np.argsort(-post)[: 4 * PAIR_TOPS] if grid.g[p] != grid.h[p]]
    pairs = list(dict.fromkeys((int(grid.g[p]), int(grid.h[p])) for p in tops))[:PAIR_TOPS]
    at = {int(g): int(i) for i, (g, e, t) in enumerate(zip(grid.g, grid.ell, grid.t,
                                                          strict=True))
          if t == 0 and e == PAIR_ELL}  # fmt: skip
    scale0 = model.scale
    for g, h in pairs:
        for x in (scale0, scale0 - 0.1, scale0 - 0.2) if fit_scale else (scale0,):
            if not SCALE_GRID[0] <= x <= SCALE_GRID[-1]:
                continue
            model.set_scale(float(x))
            trial = list(pts)
            trial[0], trial[1] = at[g], at[h]
            d = depth.copy()
            d[0] = d[1] = max(float(depth[0]) / 2, LAMBDA_RANGE[0])
            ll = -np.inf
            for _ in range(2):
                for j in (0, 1):
                    d[j], ll = _best_depth(model, trial, d, j)
            if ll > best[0]:
                best = (ll, trial, d.copy(), float(x))
    model.set_scale(best[3])
    return best[1], best[2], best[0]


def _search(
    model: Model, points: list[int], lam: list[float], fit_scale: bool
) -> tuple[list[int], list[float]]:
    """Starts for K strains by coordinate ascent: each strain's placement over the
    reference nodes (t = 0, ℓ in ``SEARCH_ELL``) and its own, with its depth maximised;
    splits of two strains' summed depth (``SPLITS``), the second on the first's placement or
    a node; then the survival scale (± 0.05, ± 0.1); a few rounds. ``model`` is left at
    the scale found."""
    grid = model.grid
    nodes = np.flatnonzero((grid.t == 0) & np.isin(grid.ell, SEARCH_ELL)).tolist()
    pts, depth = list(points), np.asarray(lam, dtype=np.float64)
    current = model.loglik(pts, depth)
    best_scale = model.scale
    for _ in range(SEARCH_ROUNDS):
        before = current
        for j in range(len(pts)):
            for p in sorted(set(nodes) | {pts[j]}):
                trial = list(pts)
                trial[j] = p
                d, ll = _best_depth(model, trial, depth, j)
                if ll > current:
                    pts, current = trial, ll
                    depth[j] = d
        # splits: strain i's depth shared with strain j, on i's placement or a node (one
        # strain absorbing another is several coordinate moves from the split)
        for i in range(len(pts)):
            for j in range(len(pts)):
                if i == j:
                    continue
                total = depth[i] + depth[j]
                scale0 = model.scale
                for share, p in itertools.product(SPLITS, sorted({pts[i], *nodes})):
                    trial, split = list(pts), depth.copy()
                    trial[j] = p
                    split[i], split[j] = total * share, total * (1 - share)
                    # a split brings a second allele: the survival scale drops with it
                    for x in (scale0, *(scale0 - SCALE_DROPS)) if fit_scale else (scale0,):
                        if not SCALE_GRID[0] <= x <= SCALE_GRID[-1]:
                            continue
                        model.set_scale(float(x))
                        ll = model.loglik(trial, split)
                        if ll > current:
                            pts, depth, current, best_scale = trial, split, ll, float(x)
                model.set_scale(best_scale)
        if model.hist is not None:  # the histogram's shape, as a plug-in
            v0 = model.hist[4]
            for v in MIX_SHAPES:
                model.set_shape(float(v))
                ll = model.loglik(pts, depth)
                if ll > current:
                    current, v0 = ll, float(v)
            model.set_shape(v0)
        if fit_scale:
            scale0 = model.scale
            for x in scale0 + np.array([-0.1, -0.05, 0.05, 0.1]):
                if not SCALE_GRID[0] <= x <= SCALE_GRID[-1]:
                    continue
                model.set_scale(float(x))
                ll = model.loglik(pts, depth)
                if ll > current:
                    current, scale0 = ll, float(x)
            model.set_scale(scale0)
            best_scale = scale0
        if current - before < 1e-6:
            break
    return pts, depth.tolist()


def _slice_scale(model: Model, c: Chain, sampler: Sampler, rng: np.random.Generator) -> None:
    """One slice-sampling update of the survival scale at β = 1 (uniform prior over
    ``SCALE_GRID``'s range); leaves ``model`` at the new scale and ``c.ll`` updated."""
    lo, hi = float(SCALE_GRID[0]), float(SCALE_GRID[-1])

    def logp(x: float) -> float:
        if not lo <= x <= hi:
            return -np.inf
        model.set_scale(x)
        return model.loglik(c.points, np.exp(c.log_lam))

    x0 = model.scale
    y = logp(x0) - rng.exponential()
    left = max(lo, x0 - 0.05 * rng.random())
    right = min(hi, left + 0.05)
    while left > lo and logp(left) > y:
        left = max(lo, left - 0.05)
    while right < hi and logp(right) > y:
        right = min(hi, right + 0.05)
    for _ in range(100):
        x1 = float(rng.uniform(left, right))
        ll = logp(x1)
        if ll > y:
            c.ll = ll
            return
        if x1 < x0:
            left = x1
        else:
            right = x1
    model.set_scale(x0)


def _rhat(x: np.ndarray) -> float:
    """Gelman-Rubin R-hat of chains x draws."""
    m, n = x.shape
    if n < 2 or m < 2:
        return float("nan")
    w = x.var(1, ddof=1).mean()
    b = n * x.mean(1).var(ddof=1)
    return float(np.sqrt(((n - 1) / n * w + b / n) / w)) if w > 0 else 1.0


def _ess(x: np.ndarray) -> float:
    """Effective sample size of chains x draws (autocorrelation summed while positive)."""
    m, n = x.shape
    if n < 4:
        return float(m * n)
    centred = x - x.mean(1, keepdims=True)
    var = (centred**2).mean()
    if var == 0:
        return float(m * n)
    rho_sum = 0.0
    for lag in range(1, n // 2):
        rho = float((centred[:, lag:] * centred[:, :-lag]).mean() / var)
        if rho <= 0:
            break
        rho_sum += rho
    return float(m * n / (1 + 2 * rho_sum))


def fit_mixture(
    sp: SpeciesPanel,
    terms: UnitTerms,
    one: OneStrain,
    k: int,
    *,
    hist: Histogram | None = None,
    others: Others | None = None,
    seed: int = 0,
    temperatures: int = TEMPERATURES,
    sweeps: int = SWEEPS,
    burn: int = BURN,
    chains: int = CHAINS,
    fixed_scale: float | None = None,
) -> MixtureFit:
    """K strains by MCMC (module docstring), started from the K = 1 fit ``one`` and the
    residuals; background and survival scale at ``one``'s."""
    rng = np.random.default_rng(seed)
    others = Others.none(len(sp.units)) if others is None else others
    scale = float(one.extra["survival_scale"] if fixed_scale is None else fixed_scale)
    nc = _not_carried(terms, others, one.beta1, int(one.extra["slab1"]))
    model = Model(sp, terms, one.grid, nc, hist, scale)
    # starts and proposals: the K = 1 fit, then each residual (a strain fitted with the
    # ones before at their posterior means)
    posts = [one.point_posterior()]
    starts, centres = [int(np.argmax(posts[0]))], [np.log(one.depth_summary()[0])]
    acc = Others(others.p.copy(), others.depth.copy(), others.f.copy())
    fitted = one
    for _ in range(1, k):
        p_new = fitted.unit_carriage
        acc = Others(
            1 - (1 - acc.p) * (1 - p_new),
            acc.depth + p_new * fitted.unit_depth,
            np.clip(1 - (1 - acc.f) * (1 - p_new * fitted.unit_f), 0.0, 1.0),
        )
        fitted = fit_one(sp, terms, acc, one.beta1, scale)
        posts.append(fitted.point_posterior())
        starts.append(int(np.argmax(posts[-1])))
        centres.append(np.log(fitted.depth_summary()[0]))
    if k > 1:  # one strain absorbs others (by raising survival): search the nodes too
        pts, depth, _ = _pair_splits(model, one, starts, np.exp(np.asarray(centres)),
                                     fixed_scale is None)  # fmt: skip
        starts, centres = _search(model, pts, list(depth), fixed_scale is None)
        centres = list(np.log(centres))
        nodes = np.zeros(len(one.grid.g))
        nodes[starts] = 1.0
        posts.append(nodes / nodes.sum())
    prior = np.exp(one.grid.log_prior - one.grid.log_prior.max())
    q = (1 - PROPOSAL_PRIOR) * np.mean(posts, 0) + PROPOSAL_PRIOR * prior / prior.sum()
    ref = Reference(q / q.sum(), np.asarray(centres))
    sampler = Sampler(model, k, ref, rng)
    chain_list = []
    for _ in range(chains):
        log_lam = np.clip(np.asarray(centres) + rng.normal(0, 0.05, k), *np.log(LAMBDA_RANGE))
        chain_list.append(Chain(list(starts), log_lam, model.loglik(starts, np.exp(log_lam))))
    if k > 1 and fixed_scale is None:
        # One strain absorbs another by raising its survival (more alleles, more k-mers):
        # the scale is refitted for K strains, as the posterior median of a pilot run at
        # β = 1 that samples it too (a plug-in, as for K = 1).
        scales = []
        for chain in chain_list:  # one shared scale, carried from chain to chain
            chain.ll = model.loglik(chain.points, np.exp(chain.log_lam))
            for i in range(PILOT):
                sampler.sweep(chain, 1.0)
                _slice_scale(model, chain, sampler, rng)
                if i >= PILOT // 2:
                    scales.append(model.scale)
        model.set_scale(float(np.median(scales)))
        for chain in chain_list:
            chain.ll = model.loglik(chain.points, np.exp(chain.log_lam))
    betas = (np.arange(temperatures + 1) / temperatures) ** 3

    def log_ratio(points: list[int], log_lam: np.ndarray, ll: float) -> float:
        """log (L prior / r): the path's increment per unit β."""
        return ll + sampler.log_prior(points, log_lam) - ref.log_density(points, log_lam)

    # rung 0: the reference, sampled exactly
    draws = [ref.draw(rng) for _ in range(chains * sweeps)]
    increments = np.array([
        log_ratio(pts, lam, model.loglik(pts, np.exp(lam)))
        if np.isfinite(_log_lambda_prior(lam).sum()) else -np.inf
        for pts, lam in draws
    ])  # fmt: skip
    log_z = float(logsumexp(betas[1] * increments) - np.log(len(increments)))
    trace_ll = np.empty((chains, sweeps))
    trace_inc = np.empty((chains, sweeps))
    trace_lam = np.empty((chains, sweeps, k))
    trace_points = np.empty((chains, sweeps, k), dtype=np.int64)
    for j in range(1, temperatures + 1):
        beta = float(betas[j])
        for c, chain in enumerate(chain_list):
            for _ in range(burn):
                sampler.sweep(chain, beta)
            for i in range(sweeps):
                sampler.sweep(chain, beta)
                order = np.argsort(chain.log_lam)
                trace_ll[c, i] = chain.ll
                trace_inc[c, i] = log_ratio(chain.points, chain.log_lam, chain.ll)
                trace_lam[c, i] = chain.log_lam[order]
                trace_points[c, i] = np.asarray(chain.points)[order]
        if j < temperatures:
            step = float(betas[j + 1] - beta)
            log_z += float(logsumexp(step * trace_inc.ravel()) - np.log(trace_inc.size))
    log_bf = log_z - one.extra["log_k0"]  # against K = 0, as the K = 1 fit
    # posterior summaries at β = 1
    flat_points = trace_points.reshape(-1, k)
    flat_lam = np.exp(trace_lam.reshape(-1, k))
    pick = np.linspace(0, len(flat_lam) - 1, min(SAMPLES, len(flat_lam))).astype(int)
    summaries = [model.carriage(list(flat_points[i]), flat_lam[i]) for i in pick]
    carriage = np.mean([c for c, _ in summaries], 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        depth = np.where(carriage > 0, np.mean([d for _, d in summaries], 0) / carriage, 0.0)
    diagnostics = {
        "survival_scale": model.scale,
        **({"histogram_shape": model.hist[4]} if model.hist is not None else {}),
        "rhat_ll": _rhat(trace_ll),
        "ess_ll": _ess(trace_ll),
        **{f"rhat_log_lambda{i + 1}": _rhat(trace_lam[:, :, i]) for i in range(k)},
        **{f"ess_log_lambda{i + 1}": _ess(trace_lam[:, :, i]) for i in range(k)},
    }
    return MixtureFit(k, log_bf, flat_points[pick], flat_lam[pick], carriage, depth,
                      diagnostics)  # fmt: skip


def strain_number_posterior(log_bf: dict[int, float], pi: float) -> dict[int, float]:
    """P(K = k | data) for k = 0 .. max: K ≥ 1 with prior ``pi``, then a truncated
    geometric (``STRAIN_RATIO``) over k; ``log_bf[k]`` against K = 0 (``log_bf[0]`` = 0)."""
    ks = sorted(set(log_bf) | {0})
    top = max(ks)
    geo = np.array([STRAIN_RATIO ** (k - 1) for k in range(1, top + 1)])
    geo /= geo.sum()
    log_p = np.array([np.log1p(-pi) if k == 0 else np.log(pi) + np.log(geo[k - 1]) + log_bf[k]
                      for k in ks])  # fmt: skip
    p = np.exp(log_p - logsumexp(log_p))
    return dict(zip(ks, p.tolist(), strict=True))
