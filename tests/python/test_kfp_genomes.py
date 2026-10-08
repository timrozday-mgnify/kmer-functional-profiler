"""kfp-genomes (phase 12): evidence layer, panel, place and update."""

import numpy as np
import polars as pl
import pytest
from scipy.stats import chi2

from kfp_genomes.evidence import (
    Evidence,
    Histogram,
    UnitTerms,
    log_likelihood,
    mixture_loglik,
    zero_hit_present,
)
from kmer_functional_profiler.query import coverage_interval, ztnb_one_rate

# --- Evidence layer (step 1) ---


def fit_terms(m: float, c: float, pi: float) -> tuple[float, float]:
    """k and H of a unit fitted at coverage_zi c and present_zi pi."""
    return pi * m * -np.expm1(-c), pi * m * c


def test_likelihood_peaks_at_coverage_zi_and_present_zi() -> None:
    m, c, pi = 200.0, 2.0, 0.6
    k, h = fit_terms(m, c, pi)
    d = np.linspace(0.5, 8, 1501)[:, None]
    f = np.linspace(0.05, 1, 951)[None, :]
    ll = log_likelihood(k, h, m, d, f)
    i, j = np.unravel_index(np.argmax(ll), ll.shape)
    assert d[i, 0] == pytest.approx(c, abs=0.01) and f[0, j] == pytest.approx(pi, abs=0.002)


def test_profile_interval_over_survival_is_coverage_interval() -> None:
    m, c, pi, w = 200.0, 3.0, 0.3, 2.5
    k, h = fit_terms(m, c, pi)
    d = np.geomspace(0.5, 20, 40001)
    f_hat = np.minimum(1.0, k / (m * -np.expm1(-d)))  # the maximum over f at each D
    drop = log_likelihood(k, h, m, c, pi, w) - log_likelihood(k, h, m, d, f_hat, w)
    inside = d[drop <= chi2.ppf(0.95, 1) / 2]
    lo, hi = coverage_interval(np.array([c]), np.array([k]), np.array([w]))
    assert inside.min() == pytest.approx(lo[0], rel=1e-3)
    assert inside.max() == pytest.approx(hi[0], rel=1e-3)


def test_zero_hit_term_matches_simulation() -> None:
    rng = np.random.default_rng(1)
    m, f, d, n = 40, 0.4, 0.1, 200_000
    present = rng.binomial(m, f, n)
    p_zero = (rng.poisson(d * present) == 0).mean()
    expect = np.exp(log_likelihood(0, 0, m, d, f))
    assert p_zero == pytest.approx(expect, abs=4 * np.sqrt(expect * (1 - expect) / n))
    assert np.exp(log_likelihood(0, 0, m, d, 1.0)) == pytest.approx(np.exp(-d * m))  # f = 1: e^-Dm
    eps = zero_hit_present(np.array([5.0, 500.0]))
    assert 1 > eps[0] > eps[1] > 0  # more k-mers, less chance to be missed when present


def test_histogram_with_one_rate_is_rate_mixtures_one_rate_fit() -> None:
    rng = np.random.default_rng(2)
    rows = []
    for unit, rate in ((7, 3.0), (9, 8.0)):
        hits = rng.poisson(rate, 400)
        hits = hits[hits > 0][:60]
        values, counts = np.unique(hits, return_counts=True)
        rows += [(unit, int(v), int(c)) for v, c in zip(values, counts, strict=True)]
    hist = Histogram.from_frame(pl.DataFrame(rows, schema=["unit", "hits", "kmers"], orient="row"))
    mu, ll = ztnb_one_rate(hist.row, hist.h, hist.c, hist.v)
    assert mixture_loglik(hist, mu[:, None]) == pytest.approx(ll)
    # and the fit's rate maximises it
    for scale in (0.9, 1.1):
        assert (mixture_loglik(hist, scale * mu[:, None]) < ll).all()


def toy_profile() -> pl.DataFrame:
    """Two hit units fitted at coverage 2 (one explained away), the profile's columns."""
    m = np.array([100.0, 80.0, 60.0])
    return pl.DataFrame({
        "unit": [3, 5, 8], "hits": [150, 90, 4], "reads": [100, 60, 4], "m_g": m,
        "coverage_zi": [2.0, 1.5, 0.0], "present_zi": [0.8, 0.9, 0.0],
        "coverage_zi_dispersion": [1.0, 1.2, 1.0], "present_prob": [0.999, 0.99, 0.0],
        "present_llr": [12.0, 10.0, None],
    })  # fmt: skip


def test_unit_ratios_favour_the_fitted_depth_and_doubt_missed_units() -> None:
    ev = Evidence.from_profile(toy_profile())
    assert ev.index(np.array([5, 4, 8])).tolist() == [1, -1, 2]
    assert ev.informative.tolist() == [True, True, False]
    terms = UnitTerms.build(ev, np.array([3, 5, 8, 11]), np.array([100.0, 80.0, 60.0, 50.0]))
    assert terms.hit.tolist() == [True, True, False, False] and terms.blank[2]
    d = np.array([[2.0], [0.2], [2.0], [2.0]])  # unit 3 at its depth and well below
    lr = terms.log_ratio(np.c_[d, d * 5], np.full((4, 2), 0.8), beta=0.01)
    assert lr[0, 0] > 0 > lr[0, 1]  # at its fitted depth, carriage beats not carried
    assert lr[2].tolist() == [0.0, 0.0]  # explained away: no evidence either way
    # a zero-hit unit: carried at depth 2 is unlikely, against no hits from elsewhere
    expect = log_likelihood(0, 0, 50, 2.0, 0.8) - np.log1p(-0.01 * (1 - terms.eps[3]))
    assert lr[3, 0] == pytest.approx(expect) and lr[3, 0] < -30
