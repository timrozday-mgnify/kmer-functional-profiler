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
"""

from dataclasses import dataclass
from functools import cached_property
from typing import Any, Final, Self

import numpy as np
from scipy.stats import gamma

T_GRID: Final = np.concatenate([[0.0], np.geomspace(1e-4, 50.0, 1199)])  # divergence
MAX_LAG: Final = 300  # co-survival is tabulated for windows up to this far apart
CATEGORIES: Final = 16


@dataclass(frozen=True)
class SurvivalModel:
    """The survival model of ``aai`` (see the module docstring) for k-mers of length ``k``."""

    k: int
    shape: float | None = None
    region: float = 1.0
    categories: int = CATEGORIES

    @classmethod
    def from_json(cls, model: dict[str, Any] | None, k: int) -> Self:
        """The model of an ``aai_model.json`` (None: independent substitutions)."""
        if model is None:
            return cls(k)
        return cls(k, model["shape"], model["region"], model.get("categories", CATEGORIES))

    @property
    def rates(self) -> np.ndarray:
        """The categories' rates: gamma quantiles at the categories' midpoints, mean 1."""
        if self.shape is None:
            return np.ones(1)
        q = (np.arange(self.categories) + 0.5) / self.categories
        r = gamma.ppf(q, self.shape, scale=1 / self.shape)
        return r / r.mean()

    @property
    def switch(self) -> float:
        """Chance per site that the rate category is redrawn."""
        return 0.0 if self.shape is None else 1 / self.region

    @cached_property
    def tables(self) -> dict[str, np.ndarray]:
        """Over ``T_GRID``: ``a``, ``S``; ``rho`` (T x ``MAX_LAG``), the correlation
        (B_j - S^2) / (S (1 - S)) of the survival of windows j = 1.. apart; and its running
        sums ``R1`` = sum_{i<=j} rho_i and ``R2`` = sum_{i<=j} i rho_i, for :meth:`overlap`."""
        k, t = self.k, T_GRID[:, None]
        r, rho = self.rates, self.switch
        pi = np.full(len(r), 1 / len(r))
        e = np.exp(-t * r)  # (T, C): a site unchanged, by category

        def right(v: np.ndarray) -> np.ndarray:  # v @ P: stay, or redraw from pi
            return (1 - rho) * v + rho * v.sum(axis=1, keepdims=True) * pi

        def left(x: np.ndarray) -> np.ndarray:  # P @ x
            return (1 - rho) * x + rho * (x @ pi)[:, None]

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
        j = np.arange(1, MAX_LAG + 1)
        return {"a": a, "S": survive, "rho": corr,
                "R1": corr.cumsum(axis=1), "R2": (j * corr).cumsum(axis=1)}  # fmt: skip

    def _at(self, x: np.ndarray, key: str, out: str) -> np.ndarray:
        """``out`` at the grid point where ``key`` (a or S, both decreasing in t) is ``x``."""
        tab = self.tables
        return np.interp(np.asarray(x, dtype=np.float64), tab[key][::-1], tab[out][::-1])

    def survival(self, a: np.ndarray) -> np.ndarray:
        """S at identity ``a``."""
        return self._at(a, "a", "S")

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
        return (1 - w) * self.overlap_at(lo, keep, windows) + w * self.overlap_at(hi, keep, windows)


def check_aai_model(model: dict[str, Any], params: dict[str, Any]) -> None:
    """Raise ValueError unless ``model`` is a usable survival model for an index built with
    ``params``: ``survival`` "markov_gamma", ``shape`` > 0, ``region`` >= 1, and
    ``k``/``alphabet``, if it records them (the aai-model workflow's fit does), equal to the
    index's."""
    if model.get("survival") != "markov_gamma":
        raise ValueError('survival must be "markov_gamma"')
    if not model.get("shape", 0) > 0 or not model.get("region", 0) >= 1:
        raise ValueError("shape must be > 0 and region >= 1")
    differ = [p for p in ("k", "alphabet") if p in model and model[p] != params.get(p)]
    if differ:
        raise ValueError(f"model fitted for another {', '.join(differ)}")
