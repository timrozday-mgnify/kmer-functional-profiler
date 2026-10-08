"""``place`` (phase 12, steps 3-5): the genome profile from the functional profile.

A sample holds K_s ≥ 0 strains of species s. Strain k has depth λ_k and a placement
θ = (g, h, t, ℓ): position t on the neighbour-graph edge from reference g to its neighbour h,
pulled a fraction ℓ towards the species average. Each unit copies from g, from h or from
the species average (weights (1 - ℓ)(1 - t), (1 - ℓ) t, ℓ: a two-reference Li-Stephens
copying model without linkage). From source σ the unit is carried with probability x̃_σ
(a reference's completeness-corrected carriage, or the species' q), at depth λ n_σ and
survival f_σ. Carriage is summed out per unit, so for one strain the likelihood is

    L(θ, λ) = Π_u Σ_σ w_σ(θ) [x̃_σu L_u(carried at λ n_σu, f_σu) + (1 - x̃_σu) L_u(not)]

with the unit likelihoods of :mod:`kfp_genomes.evidence`. Units other species carry enter
through those species at their posterior means (per-component rounds).

*One strain:* exact enumeration of placements (reference nodes; each neighbour pair at
t = ¼, ½, ¾; ℓ on a grid weighted by Beta(1, 4); ℓ = 1 is one point) and of λ on a log grid
weighted by a log-normal prior. This gives the posterior over placements and the marginal
likelihood against K = 0.

*Background:* a unit not carried by the species' strains is present from outside the panel
with prior β, at a depth log-normal around a slab centre (evidence.py). β and the centre are
fitted per species by marginal likelihood (type II ML), separately under K = 0 and K = 1, on
a grid: an organism outside the panel that shares many of a species' units, at its own
depth, is the background's, rather than being taken for the species.

*Survival:* a sample strain's alleles keep a share of the unit's k-mers that is its
references' times a survival scale (strains diverge from their references genome-wide),
fitted per species by marginal likelihood like the background. A fixed reference survival
cost a strain 10% more diverged than its references ~2 nats per unit.

*Presence:* P(K_s ≥ 1) from the Bayes factor and a prior π fitted by empirical Bayes over the
candidates (with a Beta(1, 9) hyperprior: most screened species are absent).
"""

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Final

import numpy as np
import polars as pl
from scipy.special import expit, logsumexp
from scipy.stats import norm

from kfp_genomes.evidence import BG_DEPTHS, F_FLOOR, Evidence, UnitTerms
from kfp_genomes.panel import Panel, SpeciesPanel
from kmer_functional_profiler.query import gather

T_INTERIOR: Final = np.array([0.25, 0.5, 0.75])
T_NODE_WEIGHT: Final = 1 / 8  # trapezoid on t ∈ {0, ¼, ½, ¾, 1}: ends ⅛, inside ¼
T_WEIGHT: Final = 1 / 4
ELL_GRID: Final = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 1.0])
ELL_SHAPE: Final = 4.0  # ℓ ~ Beta(1, ELL_SHAPE): most strains near a reference
LAMBDA_RANGE: Final = (0.01, 1000.0)
LAMBDA_COARSE: Final = 25  # log grid over LAMBDA_RANGE
LAMBDA_FINE: Final = 21  # points around the mode (replacing the coarse ones there)
ZOOMS: Final = 8  # refinements of λ's grid at most
LAMBDA_SD: Final = 1.5  # log-normal prior on λ: median 1x, sd 1.5 in ln
# A sample strain may carry a unit at other copies than its reference (copy-number
# variation): with probability COPY_CHANGE, uniform over COPIES.
COPY_CHANGE: Final = 0.05
COPIES: Final = (1, 2, 3, 4)
SCALE_GRID: Final = np.round(np.arange(0.5, 1.301, 0.05), 2)  # strain survival / reference's
SCALE_ROUNDS: Final = 4  # refits of the scale on the posterior's points
BETA_GRID: Final = expit(np.linspace(-9, 4.5, 28))  # 1e-4 .. 0.99, even in log-odds
PRESENT_A, PRESENT_B = 1.0, 9.0  # Beta hyperprior on π, the share of candidates present
CORE: Final = 0.9  # q from which a unit is core (the candidate screen)
MIN_CORE_UNITS: Final = 5  # screen: core units detected (present_prob >= 0.5)
MIN_CORE_FRACTION: Final = 0.1  # screen: share of core units detected
ROUNDS: Final = 3  # per-component rounds (the others at their posterior means)
NEGLIGIBLE: Final = 30.0  # log posterior below the best at which a λ point is skipped
KEEP_POINTS: Final = 1e-6  # posterior mass from which a point is kept for unit summaries
BATCH: Final = 2**23  # cells (points x units x λ) per evaluation batch
NOVEL: Final = 0.5  # ℓ above which a strain is reported novel (P(ℓ > NOVEL) > 0.5)


def ell_weights() -> np.ndarray:
    """Prior mass of each ℓ grid point: Beta(1, 4) over the cells between midpoints."""
    edges = np.r_[0.0, (ELL_GRID[1:] + ELL_GRID[:-1]) / 2, 1.0]
    cdf = 1 - (1 - edges) ** ELL_SHAPE
    return np.asarray(np.diff(cdf))


def lambda_grid(
    center: float | None = None, half_width: float = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    """λ points and their log quadrature weights under λ's prior (log-normal, truncated to
    ``LAMBDA_RANGE``; trapezoid on log λ): a log grid over the whole range, or
    ``LAMBDA_FINE`` points evenly over log ``center`` ± ``half_width`` only (the mass
    outside such a window is taken as negligible: the caller's zoom keeps ± 8 sd)."""
    lo, hi = np.log(LAMBDA_RANGE)
    if center is None or half_width <= 0:
        x = np.linspace(lo, hi, LAMBDA_COARSE)
    else:
        c = np.log(center)
        x = np.linspace(max(c - half_width, lo), min(c + half_width, hi), LAMBDA_FINE)
    width = np.full(len(x), x[1] - x[0] if len(x) > 1 else 1.0)
    width[[0, -1]] /= 2
    mass = norm.cdf(hi, 0.0, LAMBDA_SD) - norm.cdf(lo, 0.0, LAMBDA_SD)
    return np.exp(x), norm.logpdf(x, 0.0, LAMBDA_SD) + np.log(width) - np.log(mass)


def _zoom(lam: np.ndarray, log_post: np.ndarray) -> tuple[float, float, bool]:
    """Centre and half-width (on log λ) for the next grid, and whether the current grid
    resolves the posterior: then the mean ± 8 sd, else three steps either side of the
    mode."""
    p = np.exp(log_post - logsumexp(log_post))
    x = np.log(lam)
    mean = float((p * x).sum())
    sd = float(np.sqrt(max((p * (x - mean) ** 2).sum(), 0.0)))
    step = float(x[1] - x[0]) if len(x) > 1 else 1.0
    if sd < 2 * step:  # unresolved: a few points carry it
        return float(np.exp(x[int(np.argmax(log_post))])), 3 * step, False
    return float(np.exp(mean)), 8 * sd, True


@dataclass(frozen=True)
class Placements:
    """The placement grid of one species: per point (g, h, t, ℓ) and its log prior."""

    g: np.ndarray
    h: np.ndarray
    t: np.ndarray
    ell: np.ndarray
    log_prior: np.ndarray

    @classmethod
    def build(cls, sp: SpeciesPanel) -> "Placements":
        n = len(sp.genomes)
        ell_w = ell_weights()
        node = np.zeros(n)
        pair: dict[tuple[int, int], float] = {}
        for (g, h), w in zip(sp.edges.tolist(), sp.edge_prior.tolist(), strict=True):
            node[g] += w * T_NODE_WEIGHT
            node[h] += w * T_NODE_WEIGHT
            if g != h:  # the edge g -> h at t is the edge h -> g at 1 - t: one point
                key = (min(g, h), max(g, h))
                pair[key] = pair.get(key, 0.0) + w * T_WEIGHT
            else:  # a lone genome: its node takes the whole edge
                node[g] += w * (1 - 2 * T_NODE_WEIGHT)
        gs, hs, ts, ws = list(range(n)), list(range(n)), [0.0] * n, node.tolist()
        for (g, h), w in pair.items():
            gs += [g] * len(T_INTERIOR)
            hs += [h] * len(T_INTERIOR)
            ts += T_INTERIOR.tolist()
            ws += [w] * len(T_INTERIOR)
        base = np.array(ws)
        base /= base.sum()
        k = len(ELL_GRID) - 1  # ℓ < 1 for each point, then one ℓ = 1 point
        g = np.r_[np.repeat(gs, k), 0]
        h = np.r_[np.repeat(hs, k), 0]
        t = np.r_[np.repeat(ts, k), 0.0]
        ell = np.r_[np.tile(ELL_GRID[:-1], len(gs)), 1.0]
        prior = np.r_[np.repeat(base, k) * np.tile(ell_w[:-1], len(gs)), ell_w[-1]]
        return cls(g.astype(np.int64), h.astype(np.int64), t, ell, np.log(prior))

    def weights(self, at: np.ndarray | slice = slice(None)) -> tuple[np.ndarray, ...]:
        """Copy weights of g, h and the species average."""
        t, ell = self.t[at], self.ell[at]
        return (1 - ell) * (1 - t), (1 - ell) * t, ell


@dataclass(frozen=True)
class Others:
    """Per unit of a species, what the other species of its component contribute: the
    chance one carries it, their summed depth and survival when they do."""

    p: np.ndarray
    depth: np.ndarray
    f: np.ndarray

    @classmethod
    def none(cls, n: int) -> "Others":
        return cls(np.zeros(n), np.zeros(n), np.zeros(n))


@dataclass
class OneStrain:
    """A species fitted with K ≤ 1 (exact grid)."""

    species: int
    log_bf: float  # log P(data | K = 1) / P(data | K = 0)
    beta0: float
    beta1: float
    lam: np.ndarray  # λ grid
    log_post: np.ndarray  # points x λ, normalised over K = 1
    grid: Placements
    unit_carriage: np.ndarray  # U: P(carried | data, K = 1)
    unit_depth: np.ndarray  # U: E[λ n | carried]
    unit_f: np.ndarray  # U: E[f | carried]
    present_prob: float = 0.0
    extra: dict[str, float] = field(default_factory=dict)

    def lambda_posterior(self) -> np.ndarray:
        return np.asarray(np.exp(logsumexp(self.log_post, axis=0)))

    def _cells(self) -> tuple[np.ndarray, np.ndarray]:
        """λ's posterior as mass spread evenly over cells of log λ around each grid point
        (edges half-way between points): the edges and the cumulative mass at them."""
        x = np.log(self.lam)
        mid = (x[1:] + x[:-1]) / 2
        step = x[1] - x[0] if len(x) > 1 else 0.0
        edges = np.r_[x[0] - step / 2, mid, x[-1] + step / 2]
        return edges, np.r_[0.0, np.cumsum(self.lambda_posterior())]

    def depth_summary(self) -> tuple[float, float, float]:
        """Posterior mean of λ and its 95% interval (λ's cells, :meth:`_cells`)."""
        p = self.lambda_posterior()
        edges, cdf = self._cells()
        lo, hi = np.interp([0.025, 0.975], cdf / cdf[-1], edges)
        return float((p * self.lam).sum()), float(np.exp(lo)), float(np.exp(hi))

    def sample(self, n: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
        """``n`` posterior draws: placement points and λ (uniform within λ's cell)."""
        flat = np.exp(self.log_post - self.log_post.max()).ravel()
        cells = rng.choice(flat.size, n, p=flat / flat.sum())
        point, j = np.divmod(cells, self.log_post.shape[1])
        edges, _ = self._cells()
        return point, np.exp(rng.uniform(edges[j], edges[j + 1]))

    def point_posterior(self) -> np.ndarray:
        return np.asarray(np.exp(logsumexp(self.log_post, axis=1)))


def survival(sp: SpeciesPanel, scale: float) -> np.ndarray:
    """Per source (references, then the species average) and unit, the sample allele's
    survival: the source's, times the strain's survival ``scale``, at most 1."""
    return np.asarray(np.clip(np.vstack([sp.f, sp.f_mean[None]]) * scale, F_FLOOR, 1.0))


def _carried(
    sp: SpeciesPanel,
    terms: UnitTerms,
    others: Others,
    lam: np.ndarray,
    scale: float = 1.0,
    sources: slice = slice(None),
) -> np.ndarray:
    """log L(carried) per source σ (G references, then the species average; ``sources`` of
    them), unit and λ, relative to M for hit units, with the other species at their means
    and survival times ``scale``."""
    n_src = np.vstack([sp.n, sp.n_mean[None]])[sources]  # S x U
    f_src = survival(sp, scale)[sources]
    rows = np.arange(len(sp.units))
    out = np.empty((len(n_src), len(rows), len(lam)))
    lp = np.log(np.maximum(others.p, 1e-300))[:, None]
    l1p = np.log1p(-np.minimum(others.p, 1 - 1e-12))[:, None]
    change = copy_change(terms, lam, survival(sp, scale)[-1], rows)
    for s in range(len(n_src)):
        d = lam[None, :] * n_src[s][:, None]
        f = np.broadcast_to(f_src[s][:, None], d.shape)
        out[s] = np.logaddexp(np.log1p(-COPY_CHANGE) + terms.log_carried(d, f, rows), change)
        if others.p.any():  # the others' allele as well (independent k-mer survival)
            od = np.broadcast_to(others.depth[:, None], d.shape)
            of = np.broadcast_to(others.f[:, None], d.shape)
            both = terms.log_carried_alleles([d, od], [f, of], rows)
            out[s] = np.logaddexp(lp + both, l1p + out[s])
    return out


def copy_change(
    terms: UnitTerms, lam: np.ndarray, f: np.ndarray, rows: np.ndarray | slice = slice(None)
) -> np.ndarray:
    """log of ``COPY_CHANGE`` x L(carried) with the copies changed from the reference's:
    uniform over ``COPIES`` (survival ``f``, the species' mean), per unit and λ."""
    lam = np.atleast_1d(lam)
    f = np.asarray(f)[rows] if np.ndim(f) else f
    parts = [terms.log_carried(lam[None, :] * c, np.broadcast_to(np.reshape(f, (-1, 1)),
                                                                 (np.size(f), len(lam))), rows)
             for c in COPIES]  # fmt: skip
    return np.asarray(np.log(COPY_CHANGE) + logsumexp(parts, axis=0) - np.log(len(COPIES)))


def _not_carried(terms: UnitTerms, others: Others, beta: float, slab: int) -> np.ndarray:
    """log L(not carried by the species) per unit: the others carry it, or the
    background (prior ``beta``, depth around ``BG_DEPTHS[slab]``)."""
    nc = terms.log_not_carried(beta, slab)
    if not others.p.any():
        return nc
    rows = np.arange(len(nc))
    alone = terms.log_carried(others.depth[:, None], others.f[:, None], rows)[:, 0]
    lp = np.log(np.maximum(others.p, 1e-300))
    return np.asarray(np.logaddexp(lp + alone, np.log1p(-np.minimum(others.p, 1 - 1e-12)) + nc))


def _log_a(sp: SpeciesPanel, carried: np.ndarray, not_carried: np.ndarray) -> np.ndarray:
    """log a_σ = log(x̃ L(carried) + (1 - x̃) L(not)) per source, unit and λ."""
    xt = np.vstack([sp.xt, sp.q[None]])[:, :, None]
    with np.errstate(divide="ignore"):
        return np.asarray(
            np.logaddexp(
                np.log(xt) + carried, np.log1p(-np.minimum(xt, 1.0)) + not_carried[None, :, None]
            )
        )


def _evaluate(
    grid: Placements, log_a: np.ndarray, at: np.ndarray, lam_at: np.ndarray | slice = slice(None)
) -> np.ndarray:
    """Σ_u log Σ_σ w_σ a_σ for points ``at`` and λ columns ``lam_at`` (points x λ)."""
    a = log_a[:, :, lam_at]
    c = a.max(0)  # per unit and λ
    scaled = np.exp(a - c[None])
    cells = scaled.shape[1] * scaled.shape[2]
    out = np.empty((len(at), scaled.shape[2]))
    step = max(1, BATCH // max(cells, 1))
    avg = scaled[-1]
    for i in range(0, len(at), step):
        sel = at[i : i + step]
        wg, wh, wa = grid.weights(sel)
        mix = (wg[:, None, None] * scaled[grid.g[sel]] + wh[:, None, None] * scaled[grid.h[sel]]
               + wa[:, None, None] * avg[None])  # fmt: skip
        with np.errstate(divide="ignore"):
            out[i : i + step] = np.log(mix).sum(1)
    return np.asarray(out + c.sum(0)[None])


def _fit_background(
    sp: SpeciesPanel,
    carried: np.ndarray,
    not_carried: Callable[[float, int], np.ndarray],
    log_w_lam: np.ndarray,
) -> tuple[tuple[float, int], tuple[float, int]]:
    """The background (β, slab) under K = 0 and under K = 1 at the species average
    (ℓ = 1), each by maximum marginal likelihood over ``BETA_GRID`` x ``BG_DEPTHS``."""
    best0, best1 = (-np.inf, (0.0, 0)), (-np.inf, (0.0, 0))
    q = sp.q[:, None]
    for slab in range(len(BG_DEPTHS)):
        for beta in BETA_GRID:
            nc = not_carried(float(beta), slab)
            k0 = float(nc.sum())
            with np.errstate(divide="ignore"):
                avg = np.logaddexp(
                    np.log(q) + carried[-1], np.log1p(-np.minimum(q, 1)) + nc[:, None]
                )
            k1 = float(logsumexp(avg.sum(0) + log_w_lam))
            if k0 > best0[0]:
                best0 = (k0, (float(beta), slab))
            if k1 > best1[0]:
                best1 = (k1, (float(beta), slab))
    return best0[1], best1[1]


def _fit_scale(
    sp: SpeciesPanel, terms: UnitTerms, others: Others, grid: Placements, lam: np.ndarray,
    log_w: np.ndarray, nc: np.ndarray,
) -> float:  # fmt: skip
    """The strain's survival scale by maximum marginal likelihood over the reference nodes
    (t = 0, every ℓ; λ integrated): ``SCALE_GRID`` in steps of 0.1, then ± 0.05."""
    nodes = np.flatnonzero(grid.t == 0)

    def score(scale: float) -> float:
        log_a = _log_a(sp, _carried(sp, terms, others, lam, scale), nc)
        ll = _evaluate(grid, log_a, nodes) + grid.log_prior[nodes, None] + log_w[None]
        return float(logsumexp(ll))

    coarse = SCALE_GRID[::2]
    best = float(coarse[int(np.argmax([score(float(x)) for x in coarse]))])
    fine = [x for x in (best - 0.05, best, best + 0.05) if SCALE_GRID[0] <= x <= SCALE_GRID[-1]]
    return float(fine[int(np.argmax([score(x) for x in fine]))])


def _fit_grid(
    sp: SpeciesPanel, terms: UnitTerms, others: Others, grid: Placements, nc: np.ndarray,
    scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:  # fmt: skip
    """λ's grid zoomed onto the posterior until it resolves it, and log (likelihood x
    prior) per (point, λ). The first (coarse) pass evaluates every node, and the pairs at
    the λ where some node is within NEGLIGIBLE of the best; later passes only the points
    within NEGLIGIBLE of the best (the rest stay -inf: negligible mass)."""
    lam, log_w = lambda_grid()
    carried = _carried(sp, terms, others, lam, scale)
    is_node = grid.t == 0
    points = np.arange(len(grid.g))
    final = False
    for i in range(ZOOMS + 1):
        log_a = _log_a(sp, carried, nc)
        ll = np.full((len(grid.g), len(lam)), -np.inf)
        if i == 0:
            nodes, pairs = np.flatnonzero(is_node), np.flatnonzero(~is_node)
            ll[nodes] = _evaluate(grid, log_a, nodes)
            best = (ll[nodes] + grid.log_prior[nodes, None]).max(0) + log_w
            live = np.flatnonzero(best >= best.max() - NEGLIGIBLE)
            if len(pairs):
                ll[np.ix_(pairs, live)] = _evaluate(grid, log_a, pairs, live)
        else:
            ll[points] = _evaluate(grid, log_a, points)
        joint = ll + grid.log_prior[:, None] + log_w[None]
        if final or i == ZOOMS:
            break
        center, half, final = _zoom(lam, logsumexp(joint, axis=0))  # resolved: once more
        per_point = logsumexp(joint, axis=1)
        points = np.flatnonzero(per_point >= per_point.max() - NEGLIGIBLE)
        lam, log_w = lambda_grid(center, half)
        carried = _carried(sp, terms, others, lam, scale)
    return lam, log_w, carried, log_a, joint


def fit_one(
    sp: SpeciesPanel,
    terms: UnitTerms,
    others: Others | None = None,
    beta: float | None = None,
    scale: float | None = None,
) -> OneStrain:
    """K ≤ 1 for one species, exactly on the grid (module docstring). ``beta`` fixes the
    background prior (both models), ``scale`` the survival scale, instead of fitting them."""
    others = Others.none(len(sp.units)) if others is None else others
    grid = Placements.build(sp)

    def not_carried(b: float, slab: int) -> np.ndarray:
        return _not_carried(terms, others, b, slab)

    def background(carried: np.ndarray) -> tuple[tuple[float, int], tuple[float, int]]:
        if beta is None:
            return _fit_background(sp, carried, not_carried, log_w)
        bg0, bg1 = _fit_background(sp, carried, lambda _, j: not_carried(beta, j), log_w)
        return (beta, bg0[1]), (beta, bg1[1])  # fixed β, the slab still fitted

    lam, log_w = lambda_grid()
    carried = _carried(sp, terms, others, lam)
    bg0, bg1 = background(carried)
    fixed_scale = scale is not None
    if scale is None:  # survival scale given the background, then the background again
        scale = _fit_scale(sp, terms, others, grid, lam, log_w, not_carried(*bg1))
        bg0, bg1 = background(_carried(sp, terms, others, lam, scale))
    beta0, beta1 = bg0[0], bg1[0]
    k0 = float(not_carried(*bg0).sum())
    nc = not_carried(*bg1)
    lam, log_w, carried, log_a, joint = _fit_grid(sp, terms, others, grid, nc, scale)
    for _ in range(0 if fixed_scale else SCALE_ROUNDS):
        # the nodes confound survival with carriage (a lower scale excuses units a node
        # carries and the strain lacks): refit it on the posterior's own points
        log_post = joint - logsumexp(joint)
        top = np.flatnonzero(logsumexp(log_post, axis=1) > np.log(KEEP_POINTS))

        def score(
            x: float, lam: np.ndarray = lam, top: np.ndarray = top, log_w: np.ndarray = log_w
        ) -> float:
            a = _log_a(sp, _carried(sp, terms, others, lam, x), nc)
            return float(logsumexp(_evaluate(grid, a, top) + grid.log_prior[top, None]
                                   + log_w[None]))  # fmt: skip

        near = [x for x in scale + np.array([-0.1, -0.05, 0.0, 0.05, 0.1])
                if SCALE_GRID[0] - 1e-9 <= x <= SCALE_GRID[-1] + 1e-9]  # fmt: skip
        best = float(near[int(np.argmax([score(float(x)) for x in near]))])
        if abs(best - scale) < 1e-9:
            break
        scale = best
        lam, log_w, carried, log_a, joint = _fit_grid(sp, terms, others, grid, nc, scale)
    log_z1 = float(logsumexp(joint))
    log_post = joint - log_z1
    keep = np.argwhere(log_post > np.log(KEEP_POINTS) + log_post.max())
    carriage, depth, surv = _unit_summaries(sp, grid, log_a, carried, lam, keep, log_post, scale)
    return OneStrain(
        species=sp.species, log_bf=log_z1 - k0, beta0=beta0, beta1=beta1, lam=lam,
        log_post=log_post, grid=grid, unit_carriage=carriage, unit_depth=depth, unit_f=surv,
        extra={"background_depth0": float(BG_DEPTHS[bg0[1]]),
               "background_depth1": float(BG_DEPTHS[bg1[1]]), "slab1": float(bg1[1]),
               "log_k0": k0, "survival_scale": scale},
    )  # fmt: skip


def _unit_summaries(
    sp: SpeciesPanel,
    grid: Placements,
    log_a: np.ndarray,
    carried: np.ndarray,
    lam: np.ndarray,
    keep: np.ndarray,
    log_post: np.ndarray,
    scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """P(carried | data), E[depth | carried] and E[f | carried] per unit, averaged over the
    kept (point, λ) cells by their posterior."""
    xt = np.vstack([sp.xt, sp.q[None]])
    n_src = np.vstack([sp.n, sp.n_mean[None]])
    f_src = survival(sp, scale)
    avg = len(n_src) - 1
    weight = np.exp(log_post[keep[:, 0], keep[:, 1]])
    weight /= weight.sum()
    p_c, e_d, e_f = (np.zeros(len(sp.units)) for _ in range(3))
    for (i, j), w in zip(keep.tolist(), weight.tolist(), strict=True):
        wg, wh, wa = (float(x[0]) for x in grid.weights(np.array([i])))
        src = [(int(grid.g[i]), wg), (int(grid.h[i]), wh), (avg, wa)]
        top = np.max([log_a[s, :, j] for s, _ in src], axis=0)
        mix = np.zeros(len(sp.units))
        c_sum, d_sum, f_sum = (np.zeros(len(sp.units)) for _ in range(3))
        for s, ws in src:
            mix += ws * np.exp(log_a[s, :, j] - top)
            c = ws * xt[s] * np.exp(carried[s, :, j] - top)
            c_sum += c
            d_sum += c * lam[j] * n_src[s]
            f_sum += c * f_src[s]
        p_c += w * c_sum / mix
        e_d += w * d_sum / mix
        e_f += w * f_sum / mix
    with np.errstate(invalid="ignore", divide="ignore"):
        return p_c, np.where(p_c > 0, e_d / p_c, 0.0), np.where(p_c > 0, e_f / p_c, 0.0)


# --- Candidates, components and presence ---


def screen(ev: Evidence, panel: Panel) -> pl.DataFrame:
    """Candidate species: core units (q ≥ ``CORE``) detected (present_prob ≥ 0.5) in at
    least ``MIN_CORE_UNITS`` and ``MIN_CORE_FRACTION`` of them; then gathered over those
    units (G0), so a species whose detected core another explains drops out. A screen
    only: presence is the model's."""
    detected = pl.DataFrame(
        {"unit": ev.unit[ev.present_prob >= 0.5].astype(np.uint32)}, schema={"unit": pl.UInt32}
    )
    core = panel.species_units.filter(pl.col("q") >= CORE).select("species", "unit")
    counts = core.group_by("species").agg(core=pl.len())
    hit = core.join(detected, on="unit", how="semi")
    passing = (
        hit.group_by("species")
        .agg(core_hit=pl.len())
        .join(counts, on="species")
        .filter(
            pl.col("core_hit") >= MIN_CORE_UNITS,
            pl.col("core_hit") >= MIN_CORE_FRACTION * pl.col("core"),
        )
    )
    if passing.height == 0:
        return passing.select("species", "core_hit", "core")
    kept = gather(
        hit.join(passing, on="species", how="semi").select(
            unit=pl.col("species").cast(pl.UInt32), hash=pl.col("unit").cast(pl.UInt64)
        ),
        np.ones(int(passing["species"].to_numpy().max()) + 1),
    ).filter(pl.col("kmers_unique") >= MIN_CORE_UNITS)
    return passing.join(kept.select(species="unit"), on="species", how="semi").sort("species")


def components(panels: dict[int, SpeciesPanel], terms: dict[int, UnitTerms]) -> list[list[int]]:
    """Candidate species linked by hit units they share (units no one's evidence touches
    do not link them)."""
    owner: dict[int, int] = {}
    parent = {s: s for s in panels}

    def find(s: int) -> int:
        while parent[s] != s:
            parent[s] = parent[parent[s]]
            s = parent[s]
        return s

    for s, sp in panels.items():
        for u in sp.units[terms[s].hit].tolist():
            if u in owner:
                parent[find(s)] = find(owner[u])
            else:
                owner[u] = s
    groups: dict[int, list[int]] = {}
    for s in panels:
        groups.setdefault(find(s), []).append(s)
    return list(groups.values())


def _others(s: int, group: list[int], fits: dict[int, OneStrain],
            panels: dict[int, SpeciesPanel]) -> Others:  # fmt: skip
    """The other species of ``s``'s component at their posterior means, on s's units."""
    units = panels[s].units
    not_p, depth, f_num = np.ones(len(units)), np.zeros(len(units)), np.zeros(len(units))
    for o in group:
        if o == s or o not in fits:
            continue
        fit, ou = fits[o], panels[o].units
        at = np.searchsorted(ou, units)
        at_c = np.minimum(at, len(ou) - 1)
        shared = ou[at_c] == units
        p = np.where(shared, fit.present_prob * fit.unit_carriage[at_c], 0.0)
        not_p *= 1 - p
        depth += p * np.where(shared, fit.unit_depth[at_c], 0.0)
        f_num += p * np.where(shared, fit.unit_f[at_c], 0.0)
    p_any = 1 - not_p
    with np.errstate(invalid="ignore", divide="ignore"):
        return Others(
            p_any,
            np.where(p_any > 0, depth / p_any, 0.0),
            np.clip(np.where(p_any > 0, f_num / p_any, 0.0), 0.0, 1.0),
        )


def presence_prior(log_bf: np.ndarray, iterations: int = 200) -> tuple[float, np.ndarray]:
    """π by empirical Bayes over the candidates (MAP under Beta(``PRESENT_A``,
    ``PRESENT_B``)) and each candidate's P(K ≥ 1)."""
    pi = PRESENT_A / (PRESENT_A + PRESENT_B)
    post = np.zeros(len(log_bf))
    for _ in range(iterations):
        post = expit(np.log(pi) - np.log1p(-pi) + log_bf)
        new = (post.sum() + PRESENT_A - 1) / (len(log_bf) + PRESENT_A + PRESENT_B - 2)
        new = float(np.clip(new, 1e-6, 1 - 1e-6))
        if abs(new - pi) < 1e-10:
            break
        pi = new
    return pi, np.asarray(post)


@dataclass
class Placed:
    """The result of :func:`place_species`: candidates, their K ≤ 1 fits and panels."""

    candidates: pl.DataFrame
    fits: dict[int, OneStrain]
    panels: dict[int, SpeciesPanel]
    terms: dict[int, UnitTerms]
    pi: float
    report: dict[str, float]


def place_species(
    profile: pl.DataFrame,
    panel: Panel,
    *,
    rounds: int = ROUNDS,
    beta: float | None = None,
    species: list[int] | None = None,
) -> Placed:
    """Screen candidates, fit each with K ≤ 1 (exact), in per-component rounds with the
    others at their posterior means, and set presence by empirical Bayes. ``species``
    skips the screen (those species only); ``beta`` fixes the background prior."""
    ev = Evidence.from_profile(profile)
    candidates = (
        screen(ev, panel)
        if species is None
        else pl.DataFrame({"species": species}, schema={"species": pl.UInt32})
    )
    panels = {s: panel.species_panel(s) for s in candidates["species"].to_list()}
    terms = {s: UnitTerms.build(ev, sp.units, sp.m_g) for s, sp in panels.items()}
    groups = components(panels, terms)
    fits: dict[int, OneStrain] = {}
    pi = PRESENT_A / (PRESENT_A + PRESENT_B)
    for r in range(max(1, rounds)):
        for group in groups:
            for s in group:
                others = _others(s, group, fits, panels) if r > 0 and len(group) > 1 else None
                before = fits[s].present_prob if s in fits else 1.0
                fits[s] = fit_one(panels[s], terms[s], others, beta)
                fits[s].present_prob = before  # until this round's presence update
        order = sorted(fits)
        pi, post = presence_prior(np.array([fits[s].log_bf for s in order]))
        for s, p in zip(order, post.tolist(), strict=True):
            fits[s].present_prob = p
        if all(len(g) == 1 for g in groups):
            break
    report = {"candidates": float(len(fits)), "components": float(len(groups)), "pi": pi}
    return Placed(candidates, fits, panels, terms, pi, report)
