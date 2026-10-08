"""Simulated profiles for kfp-genomes' tests: strains' units drawn at depth, each kept k-mer
present with the allele's survival and hit ~ Poisson, fitted as the zero-truncated EM does."""

from dataclasses import dataclass

import numpy as np
import polars as pl
from scipy.special import expit

from kmer_functional_profiler.query import ztp_lambda

M: int = 100  # kept k-mers per unit
ODDS: float = 1e-3  # the profile's prior odds of presence


@dataclass
class Strain:
    units: set[int]
    depth: float
    f: float = 0.8  # survival of every unit's allele
    copies: int = 1
    allele: int | None = None  # strains with the same allele keep the same k-mers


def _present(strain: Strain, unit: int, m: int, rng: np.random.Generator) -> np.ndarray:
    if strain.allele is None:
        return rng.random(m) < strain.f
    return np.random.default_rng([strain.allele, unit]).random(m) < strain.f


def simulate(strains: list[Strain], seed: int = 0, m: int = M) -> pl.DataFrame:
    """A profile (the columns kfp-genomes reads) of the units the strains carry."""
    return simulate_with_histogram(strains, seed, m)[0]


def simulate_with_histogram(
    strains: list[Strain], seed: int = 0, m: int = M
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """The profile and the own-k-mer histogram sidecar (every k-mer is its unit's own)."""
    rng = np.random.default_rng(seed)
    units = sorted(set().union(*(s.units for s in strains)))
    rows, hist = [], []
    for u in units:
        rate = np.zeros(m)
        for s in strains:
            if u in s.units:
                rate += s.depth * s.copies * _present(s, u, m, rng)
        hits = rng.poisson(rate)
        if (hits > 0).sum() >= 10:
            values, counts = np.unique(hits[hits > 0], return_counts=True)
            hist += [(u, int(v), int(c)) for v, c in zip(values, counts, strict=True)]
        k, h = int((hits > 0).sum()), int(hits.sum())
        if k == 0:
            continue
        c = float(ztp_lambda(np.array([h / k]))[0])
        pi = k / (m * -np.expm1(-c)) if c > 0 else 1.0
        if pi > 1 or c <= 0:  # the zero-inflated EM caps present at 1
            pi, c = 1.0, -np.log1p(-k / m) if k < m else 10.0
        llr = min(3.0 * k, 40.0) - 2.0  # a few own k-mers: some doubt
        rows.append((u, h, h, m, c, pi, 1.0, float(expit(llr + np.log(ODDS))), llr))
    profile = pl.DataFrame(
        rows,
        schema=["unit", "hits", "reads", "m_g", "coverage_zi", "present_zi",
                "coverage_zi_dispersion", "present_prob", "present_llr"],
        orient="row",
    )  # fmt: skip
    return profile, pl.DataFrame(hist, schema=["unit", "hits", "kmers"], orient="row")
