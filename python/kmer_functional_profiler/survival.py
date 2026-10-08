"""Exact k-mer survival under a Markov-modulated rate model (plan, phase 7, step 25).

A homolog differs from a protein by substitutions. Site i changes with probability
1 - exp(-r_i t) at divergence t, and a k-mer window survives (its exact k-mer is shared)
when none of its k sites changed. Rates come from ``categories`` equiprobable categories of
a gamma distribution of shape ``shape`` and mean 1, and change along the sequence as a Markov
chain: at each site the category is redrawn from the stationary distribution with
probability 1 / ``region``, the mean length of a region evolving at one rate. Then, with
E = diag(exp(-r t)) and P the chain's transition matrix:

- identity a(t) = sum_c pi_c exp(-r_c t);
- window survival S(t) = pi (E P)^(k-1) E 1;
- co-survival B_j(t) of two windows j sites apart: for j < k the survival of one window
  of k + j sites; for j >= k, (pi (E P)^(k-1) E) P^(j-k+1) beta, with beta = E (P E)^(k-1) 1
  the survival of a window from its first site's category.

iid site rates (``region`` 1) give S = a^k whatever the shape: survival exceeds a^k only
when conserved sites cluster (step 24). ``shape`` None is one category (independent
substitutions): S = a^k and B_j = a^(k + min(j, k)).

Everything depends on (``shape``, ``region``, k, t) only, so it is tabulated once on a grid
of t, and per-unit values are interpolated: ~0.06 s per model, ~0.1 s per 10^6 units.

Markov-beta (``concentration`` set; plan, phase 7, step 33): the categories are sites'
chances to stay unchanged rather than rates: the ``categories`` equal-probability slices of
Beta(a phi, (1 - a) phi), phi = ``concentration``, each at its slice's mean (so they average
exactly a), at a = exp(-t) on the same grid. A region's identity scatters around the pair's
own identity with a fixed concentration, as protein pairs show; under gamma rates the
contrast changes with t. The transfer products are the same.

``ends`` < 1 is end loss: a share 1 - ``ends`` of a unit's windows lies outside the region a
homolog aligns to and is never shared, so survival is ``ends`` S. In the correlation it
thins windows independently: rho' = ``ends`` rho (1 - S) / (1 - ``ends`` S).
ponytail: real end loss is contiguous (two ends) and varies between units, variance this
leaves out; model the aligned share per unit if intervals under-cover near identity 1.

Union (``union`` g0, ``union_slope`` g1; step 34): a unit of n members is matched against
their union. Measured at identity to the nearest member, the union still keeps more than one
member does, as m effective independent members: S_n = 1 - (1 - S)^m with
m = n^(g0 + g1 (a - ``UNION_PIVOT``)), S with end loss (the union covers members' ends too).
g0 = g1 = 0 (or n = 1) is one member. Survival depends on n, so it is per unit
(:meth:`union_at`); the tables stay one member's.

Union scatter (``spread`` c, ``spread_power`` b; step 34): a union's survival, as the query
measures it (hit k-mers over an average member's), scatters around S_n beyond window
sampling, by which members sit near the gene and how its length compares with theirs. Its
relative SD is c (1 - S_n)^b (1 - 1/n) (:meth:`union_scatter_at`): 0 for one member,
largest at low identity. Fitted on aai-model's multi-member pairs; without it, intervals of
multi-member units with many hit k-mers were far too narrow.
"""

from dataclasses import dataclass
from functools import cached_property
from typing import Any, Final, Self

import numpy as np
from scipy.stats import beta as beta_dist
from scipy.stats import gamma

T_GRID: Final = np.concatenate([[0.0], np.geomspace(1e-4, 50.0, 1199)])  # divergence
MAX_LAG: Final = 300  # co-survival is tabulated for windows up to this far apart
CATEGORIES: Final = 16
UNION_PIVOT: Final = 0.8  # identity at which the union's exponent is g0


@dataclass(frozen=True)
class SurvivalModel:
    """The survival model of ``aai`` (see the module docstring) for k-mers of length ``k``."""

    k: int
    shape: float | None = None
    region: float = 1.0
    categories: int = CATEGORIES
    concentration: float | None = None  # Markov-beta's phi, in place of ``shape``
    ends: float = 1.0
    union: float = 0.0
    union_slope: float = 0.0
    spread: float = 0.0
    spread_power: float = 0.0

    @classmethod
    def from_json(cls, model: dict[str, Any] | None, k: int) -> Self:
        """The model of an ``aai_model.json`` (None: independent substitutions)."""
        if model is None:
            return cls(k)
        extra = {
            "categories": model.get("categories", CATEGORIES),
            "ends": model.get("ends", 1.0),
            "union": model.get("union", 0.0),
            "union_slope": model.get("union_slope", 0.0),
            "spread": model.get("spread", 0.0),
            "spread_power": model.get("spread_power", 0.0),
        }
        if model["survival"] == "markov_beta":
            return cls(k, None, model["region"], concentration=model["concentration"], **extra)
        return cls(k, model["shape"], model["region"], **extra)

    @property
    def rates(self) -> np.ndarray:
        """The categories' rates: gamma quantiles at the categories' midpoints, mean 1."""
        if self.shape is None:
            return np.ones(1)
        q = (np.arange(self.categories) + 0.5) / self.categories
        r = gamma.ppf(q, self.shape, scale=1 / self.shape)
        return np.asarray(r / r.mean())

    @property
    def switch(self) -> float:
        """Chance per site that the rate category is redrawn."""
        return 0.0 if self.shape is None and self.concentration is None else 1 / self.region

    def keep_at(self, t: np.ndarray) -> np.ndarray:
        """(len(``t``) x categories): each category's chance that a site is unchanged at
        divergence ``t``."""
        t = np.asarray(t, dtype=np.float64)[:, None]
        if self.concentration is None:
            return np.asarray(np.exp(-t * self.rates))
        phi, c = self.concentration, self.categories
        a, change = np.exp(-t), -np.expm1(-t)

        def slice_means(mean: np.ndarray) -> np.ndarray:
            """The ``c`` slices' means of Beta(mean phi, (1 - mean) phi), ascending."""
            p, q = np.maximum(mean * phi, 1e-300), np.maximum((1 - mean) * phi, 1e-300)
            with np.errstate(all="ignore"):
                edges = beta_dist.ppf(np.arange(c + 1) / c, p, q)
                mass = np.diff(beta_dist.cdf(edges, p + 1, q), axis=1)  # E[X; slice] / mean
            return np.asarray(c * mean * np.nan_to_num(mass))

        # the side away from 1 keeps float precision: X for a <= 1/2, 1 - X above
        keep = np.where(a <= 0.5, slice_means(a), 1 - slice_means(change)[:, ::-1])
        return np.asarray(np.clip(keep, 0.0, 1.0))

    @cached_property
    def tables(self) -> dict[str, np.ndarray]:
        """Over ``T_GRID``: ``a``, ``S``; ``rho`` (T x ``MAX_LAG``), the correlation
        (B_j - S^2) / (S (1 - S)) of the survival of windows j = 1.. apart; and its running
        sums ``R1`` = sum_{i<=j} rho_i and ``R2`` = sum_{i<=j} i rho_i, for :meth:`overlap`."""
        k, rho = self.k, self.switch
        e = self.keep_at(T_GRID)  # (T, C): a site unchanged, by category
        pi = np.full(e.shape[1], 1 / e.shape[1])

        def right(v: np.ndarray) -> np.ndarray:  # v @ P: stay, or redraw from pi
            return np.asarray((1 - rho) * v + rho * v.sum(axis=1, keepdims=True) * pi)

        def left(x: np.ndarray) -> np.ndarray:  # P @ x
            return np.asarray((1 - rho) * x + rho * (x @ pi)[:, None])

        fwd, run = pi * e, []  # joint mass of a surviving run of sites, by the last's category
        for _ in range(2 * k - 1):
            run.append(fwd.sum(axis=1))  # runs of 1 .. 2k - 1 sites
            fwd = right(fwd) * e
        survive = run[k - 1]
        end = pi * e
        for _ in range(k - 1):
            end = right(end) * e  # a surviving window, by its last site's category
        beta = e.copy()
        for _ in range(k - 1):
            beta = e * left(beta)
        both = [run[k - 1 + j] for j in range(1, k)]  # overlapping windows
        u = end
        for _ in range(k, MAX_LAG + 1):
            u = right(u)
            both.append((u * beta).sum(axis=1))
        both_arr = np.stack(both, axis=1)[:, :MAX_LAG]
        var = survive * (1 - survive)
        cov = both_arr - survive[:, None] ** 2
        with np.errstate(invalid="ignore", divide="ignore"):
            corr = np.where(var[:, None] > 1e-12, cov / var[:, None], 0.0)
        a = e @ pi
        if self.ends < 1:  # end loss, windows thinned independently (module docstring)
            corr = corr * self.ends * (1 - survive[:, None]) / (1 - self.ends * survive[:, None])
            survive = self.ends * survive
        j = np.arange(1, MAX_LAG + 1)
        return {"a": a, "S": survive, "rho": corr,
                "R1": corr.cumsum(axis=1), "R2": (j * corr).cumsum(axis=1)}  # fmt: skip

    def _at(self, x: np.ndarray, key: str, out: str) -> np.ndarray:
        """``out`` at the grid point where ``key`` (a or S, both decreasing in t) is ``x``."""
        tab = self.tables
        return np.asarray(
            np.interp(np.asarray(x, dtype=np.float64), tab[key][::-1], tab[out][::-1])
        )

    def survival(self, a: np.ndarray) -> np.ndarray:
        """S at identity ``a`` (one member)."""
        return self._at(a, "a", "S")

    def members(self, a: np.ndarray, n_members: np.ndarray | float) -> np.ndarray:
        """The union's effective members m = n^(g0 + g1 (a - ``UNION_PIVOT``)), exponent
        floored at 0, at identity ``a`` for units of ``n_members``."""
        power = np.maximum(self.union + self.union_slope * (np.asarray(a) - UNION_PIVOT), 0.0)
        return np.asarray(np.maximum(np.asarray(n_members, dtype=np.float64), 1.0) ** power)

    def union_at(self, row: np.ndarray, n_members: np.ndarray | float) -> np.ndarray:
        """Survival of a union of ``n_members`` at integer rows ``row`` of ``T_GRID``:
        1 - (1 - S)^m. Falls with the row (as a does) when ``union_slope`` >= 0."""
        tab = self.tables
        if self.union == 0 and self.union_slope == 0:
            return np.asarray(tab["S"][row])
        m = self.members(tab["a"][row], n_members)
        return np.asarray(-np.expm1(m * np.log1p(-np.minimum(tab["S"][row], 1 - 1e-15))))

    def union_scatter_at(self, row: np.ndarray, n_members: np.ndarray | float) -> np.ndarray:
        """The relative SD of a union's survival at integer rows ``row`` of ``T_GRID``:
        ``spread`` (1 - S_n)^``spread_power`` (1 - 1/n); 0 without ``spread``."""
        n = np.maximum(np.asarray(n_members, dtype=np.float64), 1.0)
        if self.spread == 0:
            return np.zeros(np.broadcast(np.asarray(row), n).shape)
        miss = np.maximum(1 - self.union_at(row, n), 0.0)
        return np.asarray(self.spread * miss**self.spread_power * (1 - 1 / n))

    def union_survival(self, a: np.ndarray, n_members: np.ndarray | float) -> np.ndarray:
        """S_n at identity ``a`` for units of ``n_members``."""
        s = np.minimum(self.survival(a), 1 - 1e-15)
        return np.asarray(-np.expm1(self.members(a, n_members) * np.log1p(-s)))

    def identity(self, s: np.ndarray) -> np.ndarray:
        """The identity whose survival is ``s`` (beyond the grid: its end)."""
        return self._at(s, "S", "a")

    def row(self, s: np.ndarray) -> np.ndarray:
        """The fractional row of ``T_GRID`` where survival is ``s``."""
        tab = self.tables
        rows = np.arange(len(T_GRID), dtype=np.float64)
        return np.interp(np.asarray(s, dtype=np.float64), tab["S"][::-1], rows[::-1])

    def overlap_at(self, row: np.ndarray, keep: np.ndarray, windows: np.ndarray) -> np.ndarray:
        """Variance inflation d of the share of ``windows`` kept windows (sampled at rate
        ``keep``) that survive, at integer rows ``row`` of ``T_GRID``, against independent
        windows: 1 + 2 keep sum_{j < n} (1 - j / n) rho_j, n = ``windows`` (lags capped at
        ``MAX_LAG``). O(1) per unit, from the running sums."""
        tab = self.tables
        n = np.maximum(np.asarray(windows, dtype=np.float64), 1.0)
        last = np.clip(np.ceil(n).astype(int) - 2, -1, MAX_LAG - 1)  # j <= n - 1
        at = np.maximum(last, 0)
        total = np.where(last >= 0, tab["R1"][row, at] - tab["R2"][row, at] / n, 0.0)
        return 1 + 2 * np.asarray(keep, dtype=np.float64) * total

    def overlap(self, s: np.ndarray, keep: np.ndarray, windows: np.ndarray) -> np.ndarray:
        """:meth:`overlap_at` at survival ``s``, linear between grid rows."""
        s = np.atleast_1d(np.asarray(s, dtype=np.float64))
        pos = self.row(s)
        lo = np.floor(pos).astype(int)
        hi = np.minimum(lo + 1, len(T_GRID) - 1)
        keep, windows = np.broadcast_to(keep, s.shape), np.broadcast_to(windows, s.shape)
        w = pos - lo
        return np.asarray(
            (1 - w) * self.overlap_at(lo, keep, windows) + w * self.overlap_at(hi, keep, windows)
        )


def check_aai_model(model: dict[str, Any], params: dict[str, Any]) -> None:
    """Raise ValueError unless ``model`` is a usable survival model for an index built with
    ``params``: ``survival`` "markov_gamma" (``shape`` > 0) or "markov_beta"
    (``concentration`` > 0), ``region`` >= 1, ``ends`` (if given) in (0, 1], ``union`` and
    ``union_slope`` (if given) >= 0 (so survival falls with identity in every unit), and
    ``k``/``alphabet``, if it records them (the aai-model workflow's fit does), equal to the
    index's."""
    spread = {"markov_gamma": "shape", "markov_beta": "concentration"}.get(
        str(model.get("survival"))
    )
    if spread is None:
        raise ValueError('survival must be "markov_gamma" or "markov_beta"')
    if not model.get(spread, 0) > 0 or not model.get("region", 0) >= 1:
        raise ValueError(f"{spread} must be > 0 and region >= 1")
    if not 0 < model.get("ends", 1.0) <= 1:
        raise ValueError("ends must be in (0, 1]")
    if not (model.get("union", 0.0) >= 0 and model.get("union_slope", 0.0) >= 0):
        raise ValueError("union and union_slope must be >= 0")
    if not (model.get("spread", 0.0) >= 0 and model.get("spread_power", 0.0) >= 0):
        raise ValueError("spread and spread_power must be >= 0")
    differ = [p for p in ("k", "alphabet") if p in model and model[p] != params.get(p)]
    if differ:
        raise ValueError(f"model fitted for another {', '.join(differ)}")
