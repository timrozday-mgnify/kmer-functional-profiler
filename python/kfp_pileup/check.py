"""Step 6's ``aai`` check on a real sample (plan, phase 13): the profile's ``aai`` (``ani``)
against the identity of the unit's reads to its sequences, from ``units.tsv``.

No truth here: two estimators are compared, and DIAMOND's sensitivity floor and reference
bias are in the alignment side (the spike benchmark has the truth). Per breakdown (alignment
identity bin, depth tertile by ``coverage_zi``, protein length by ``len_mean``, cluster size
by ``n_members``; a breakdown whose column ``units.tsv`` lacks is skipped): units, the
median gap (``aai`` minus alignment identity) and the share of alignment identities inside
``aai_lo``-``aai_hi``. Step 6's rule flags a bin with >= ``MIN_UNITS`` units whose |median
gap| > ``MAX_GAP``, or a cover below ``MIN_COVER``: a phase-7 ``aai`` finding, not a silent
fix.
"""

from pathlib import Path
from typing import Final

import polars as pl

IDENTITY_BINS: Final = [0.8, 0.85, 0.9, 0.95, 0.97, 0.99]
LENGTH_BINS: Final = [150, 300, 500]
MEMBER_BINS: Final = [1.5, 3.5, 10.5, 100.5]
MIN_READS: Final = 10.0  # units with fewer assigned reads have too noisy an identity
MIN_UNITS: Final = 100
MAX_GAP: Final = 0.01
MIN_COVER: Final = 0.8


def aai_check(units_tsv: str | Path, min_reads: float = MIN_READS) -> pl.DataFrame:
    units = pl.read_csv(units_tsv, separator="\t", infer_schema_length=None)
    sim = "ani" if "ani" in units.columns else "aai"
    u = units.filter(
        pl.col("reads") >= min_reads, pl.col(sim).is_not_null(), pl.col("identity").is_not_null()
    )
    bins = {
        "identity_bin": pl.col("identity").cut(IDENTITY_BINS, left_closed=True),
        "depth_bin": pl.col("coverage_zi").qcut(3, allow_duplicates=True),
        "length_bin": pl.col("len_mean").cut(LENGTH_BINS),
        "members_bin": pl.col("n_members").cut(MEMBER_BINS),
    }
    sources = {"depth_bin": "coverage_zi", "length_bin": "len_mean", "members_bin": "n_members"}
    bins = {b: e for b, e in bins.items() if sources.get(b, "identity") in u.columns}
    u = u.with_columns(
        gap=pl.col(sim) - pl.col("identity"),
        inside=pl.col("identity").is_between(pl.col(f"{sim}_lo"), pl.col(f"{sim}_hi")),
        **bins,
    )
    parts = [
        u.group_by(bin=pl.col(b).cast(pl.String))
        .agg(
            units=pl.len(),
            median_gap=pl.col("gap").median(),
            cover=pl.col("inside").mean(),
            median_identity=pl.col("identity").median(),
        )
        .with_columns(breakdown=pl.lit(b.removesuffix("_bin")))
        for b in bins
    ]
    overall = u.select(
        bin=pl.lit("all"),
        units=pl.len(),
        median_gap=pl.col("gap").median(),
        cover=pl.col("inside").mean(),
        median_identity=pl.col("identity").median(),
        breakdown=pl.lit("all"),
    )
    return (
        pl.concat([overall, *parts], how="diagonal_relaxed")
        .with_columns(
            flagged=(pl.col("units") >= MIN_UNITS)
            & ((pl.col("median_gap").abs() > MAX_GAP) | (pl.col("cover") < MIN_COVER))
        )
        .select("breakdown", "bin", "units", "median_identity", "median_gap", "cover", "flagged")
        .sort("breakdown", "bin")
    )
