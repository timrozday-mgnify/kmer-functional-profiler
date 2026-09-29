"""Full-scale cost model (phase 6): index size predicted from per-cluster statistics.

``cluster_stats`` runs pass 1 of the build (distinct k-mers per unit) on a members table
and samples (hash, cluster) pairs at a low rate; it scales to the whole release when run
per cluster bucket. ``predict_cost`` applies the build's sampling rules to those counts to
give the expected postings and table sizes for any parameters, without building.
"""

from pathlib import Path
from typing import Final

import polars as pl

from kmer_functional_profiler import _core
from kmer_functional_profiler.index import (
    KEYS_PER_BUCKET_BITS,
    PIN_BITS,
    IndexParams,
    _batches,
    _kmers,
    _n_kmers,
    _smallest,
    _t_g,
    load_members,
    mask_adapters,
)

STATS_RATE: Final = 1 / 1000  # rate of the sampled (hash, cluster) pairs


def cluster_stats(
    members_path: str | Path, params: IndexParams, rate: float = STATS_RATE
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Per-cluster counts and (``hash``, ``cluster_rep``) pairs with hash <= ``rate``.

    Counts: ``cluster_rep``, ``n_members``, ``n_full_length``, ``sum_len``, ``max_len``
    and ``n_kmers`` (distinct k-mers after adapter masking, as in the build).
    """
    members = load_members(members_path)
    clusters = members.group_by("unit").agg(
        pl.col("cluster_rep").first(),
        n_members=pl.len().cast(pl.UInt32),
        n_full_length=pl.col("full_length").sum().cast(pl.UInt32),
        sum_len=pl.col("sequence").str.len_bytes().cast(pl.UInt64).sum(),
        max_len=pl.col("sequence").str.len_bytes().max().cast(pl.UInt32),
    )
    if params.mask_adapters:
        members = members.with_columns(mask_adapters(pl.col("sequence")))
    batches = list(_batches(members, params.batch_residues))
    clusters = (
        clusters.join(_n_kmers(batches, params), on="unit", how="left")
        .with_columns(pl.col("n_kmers").fill_null(0).cast(pl.UInt32))
        .sort("unit")
    )
    pairs = pl.concat(
        [_kmers(b, params, _core.max_hash(rate)).select("unit", "hash").unique() for b in batches]
    ).join(clusters.select("unit", "cluster_rep"), on="unit")
    return clusters.drop("unit"), pairs.select("hash", "cluster_rep")


def _width(n: int) -> int:
    return int(_smallest(n).itemsize)


def packed_bytes(keys: float, sets: float, values: float, max_value: int, fp_bits: int) -> float:
    """Bytes of a ``PackedTable`` with these counts, using the dtypes ``PackedTable.build``
    picks."""
    keys, sets, values = (max(int(round(n)), 1) for n in (keys, sets, values))
    buckets = 2 ** max((keys - 1).bit_length() - KEYS_PER_BUCKET_BITS, 1)
    return float(
        (buckets + 1) * _width(keys)
        + keys * (_width(2**fp_bits - 1) + _width(sets))
        + (sets + 1) * _width(values)
        + values * _width(max_value)
    )


def predict_cost(
    clusters: pl.LazyFrame, params: IndexParams, keep: float = 1.0, sets: float = 1.0
) -> dict[str, float]:
    """Expected index size from ``n_members`` and ``n_kmers`` per cluster.

    A unit expects ``t_g * n_kmers * keep`` candidates, ``keep`` being the fraction that are
    not promiscuous; floored units keep at most ``n_min`` of them. Table bytes assume one key
    per posting (an upper bound: shared k-mers share keys) and ``sets`` distinct one-value
    sets per key; the build stores each (unit, ``p_in`` level) set once, so ``sets`` < 1 is
    measured on built subsets (``tier2_sets / postings``).
    """
    candidates = pl.col("t_g") * pl.col("n_kmers") * keep
    floored = pl.col("t_g") > params.t_base
    dense_rate = pl.max_horizontal(pl.lit(params.t_dense), pl.col("t_g"))
    row = (
        clusters.with_columns(t_g=_t_g(params))
        .with_columns(
            m=pl.when(floored)
            .then(pl.min_horizontal(candidates, params.n_min))
            .otherwise(candidates)
        )
        .select(
            units=pl.len(),
            singletons=(pl.col("n_members") == 1).sum(),
            floored=floored.sum(),
            below_floor=(floored & (candidates < params.n_min)).sum(),
            kmers=pl.col("n_kmers").cast(pl.Float64).sum(),
            t_max=pl.col("t_g").max(),
            postings=pl.col("m").sum(),
            # only units with a posting are indexed, P = 1 - exp(-m) under Poisson
            with_postings=(1 - (-pl.col("m")).exp()).sum(),
            dense=(dense_rate * pl.col("n_kmers") * keep).sum()
            if params.t_dense > 0
            else pl.lit(0.0),
        )
        .collect(engine="streaming")
        .row(0, named=True)
    )
    max_value = (int(row["with_postings"]) << PIN_BITS) | (2**PIN_BITS - 1)
    n2, n_dense = row["postings"], row["dense"]
    tier2 = packed_bytes(n2, sets * n2, sets * n2, max_value, params.fp_bits)
    dense = (
        packed_bytes(n_dense, sets * n_dense, sets * n_dense, max_value, params.fp_bits)
        if params.t_dense > 0
        else 0.0
    )
    return {
        **row,
        "tier2_bytes": tier2,
        "dense_bytes": dense,
        "tier2_bytes_per_posting": tier2 / max(row["postings"], 1),
    }
