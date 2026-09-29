"""Partitioned index build (phase 6): the build of ``index.build_index`` split into jobs.

Members come in buckets that each hold whole clusters (``workflows/mgnify-subset``'s MERGE).
Everything a unit needs is local to its bucket except ``n_groups``, the number of units
holding a candidate k-mer, which feeds the score and the promiscuity cut. That is counted
exactly by a hash-partitioned reduce, with a Bloom filter of all candidate hashes to keep
the rows each bucket emits small. Units are keyed by ``cluster_rep`` until ``pack``.

1. ``candidates`` (per bucket): pass 1, the unit table and the bucket's candidate hashes.
2. ``bloom`` (once): a Bloom filter of every bucket's candidate hashes, and ``t_max``.
3. ``presence`` (per bucket): (hash, unit) of every k-mer at ``t_max`` that passes the
   filter, with its counting members and whether it is the unit's candidate, split into
   hash ranges.
4. ``groups`` (per hash range): ``n_groups`` per hash; keeps candidate rows only (which
   drops the filter's false positives) and splits them by bucket.
5. ``postings`` (per bucket): score, promiscuity cut and floor, as ``build_index``;
   ``len_cv`` and Pfam labels.
6. ``pack`` (once): all buckets' units and postings -> the index, as ``build_index``.

The result equals ``build_index`` on the concatenated members. No dense tier yet.
"""

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl
from numpy.typing import NDArray

from kmer_functional_profiler.index import (
    IndexParams,
    _kmers,
    _len_cv,
    _n_kmers,
    _prepare,
    _select_postings,
    _unit_pfam,
    _unit_table,
    finish_index,
)

BLOOM_PROBES: Final = 7
BLOOM_CHUNK: Final = 1 << 22  # hashes per vectorised probe batch
MIX: Final = np.uint64(0x9E3779B97F4A7C15)  # odd, so hash * MIX is a bijection
LOW32: Final = np.uint64(0xFFFFFFFF)
QUANTILES: Final = 1024  # candidate-hash quantiles recorded for balanced hash ranges
SAMPLE: Final = 1 << 20  # candidate hashes sampled for them
UNIT_COLUMNS: Final = ("cluster_rep", "n_members", "n_counting", "n_kmers", "t_g", "max_hash_g")


def _check(params: IndexParams) -> None:
    if params.t_dense > 0:
        raise ValueError("the partitioned build has no dense tier yet (t_dense must be 0)")


def _positions(hashes: NDArray[np.uint64], n_bits: int) -> NDArray[np.uint64]:
    """Bloom bit positions, (len(hashes), BLOOM_PROBES), by double hashing."""
    mixed = hashes * MIX  # uint64 arrays wrap
    h1, h2 = mixed >> np.uint64(32), (mixed & LOW32) | np.uint64(1)
    probes = np.arange(BLOOM_PROBES, dtype=np.uint64)
    return (h1[:, None] + probes[None, :] * h2[:, None]) % np.uint64(n_bits)


def bloom_contains(bits: NDArray[np.uint8], hashes: NDArray[np.uint64]) -> NDArray[np.bool_]:
    """Whether each hash may be in the filter (no false negatives)."""
    hashes = np.asarray(hashes, dtype=np.uint64)
    out = np.empty(len(hashes), dtype=bool)
    for i in range(0, len(hashes), BLOOM_CHUNK):
        pos = _positions(hashes[i : i + BLOOM_CHUNK], len(bits) * 8)
        out[i : i + BLOOM_CHUNK] = ((bits[pos >> np.uint64(3)] >> (pos & np.uint64(7))) & 1).all(1)
    return out


def _swap(table: pl.DataFrame, ids: pl.DataFrame, old: str, new: str) -> pl.DataFrame:
    """Replace key column ``old`` with ``new`` in place, through ``ids`` (both columns)."""
    order = [new if c == old else c for c in table.columns]
    return table.join(ids, on=old).select(order)


def _counting(members: pl.DataFrame, units: pl.DataFrame) -> pl.DataFrame:
    """Local ``unit`` ids of a bucket joined to its stage-1 unit table."""
    ids = members.group_by("unit").agg(pl.col("cluster_rep").first())
    return units.join(ids, on="cluster_rep").sort("unit")


def candidates(members_path: str | Path, prefix: str, params: IndexParams) -> None:
    """Stage 1: ``{prefix}.units.parquet``, ``.candidates.npy`` and ``.stats.json``."""
    _check(params)
    members, batches, stats = _prepare(members_path, params)
    units = _unit_table(members, _n_kmers(batches, params), params)
    thresholds = units.select("unit", "max_hash_g")
    t_max_hash = int(units["max_hash_g"].max() or 0)  # type: ignore[arg-type]
    hashes = pl.concat(
        [
            _kmers(b, params, t_max_hash)
            .join(thresholds, on="unit")
            .filter(pl.col("hash") <= pl.col("max_hash_g"))
            .select("hash")
            .unique()
            for b in batches
        ]
    ).unique()
    units.select(UNIT_COLUMNS).write_parquet(f"{prefix}.units.parquet")
    np.save(f"{prefix}.candidates.npy", np.sort(hashes["hash"].to_numpy()))
    Path(f"{prefix}.stats.json").write_text(json.dumps(stats) + "\n")


def bloom(prefixes: Sequence[str], out: str | Path, bits_per_key: float = 10.0) -> None:
    """Stage 2: ``{out}.npy`` (filter bits) and ``{out}.json`` (``t_max_hash``, sizes and
    ``QUANTILES`` + 1 quantiles of the candidate hashes, where ``presence`` cuts ranges).

    Keys are counted per bucket, so hashes shared by buckets are counted more than once and
    the filter is at most that much larger than needed.
    """
    n_keys = sum(len(np.load(f"{p}.candidates.npy", mmap_mode="r")) for p in prefixes)
    bits = np.zeros(max(int(n_keys * bits_per_key) // 8, 1), dtype=np.uint8)
    stride, sample = max(n_keys // SAMPLE, 1), []
    for p in prefixes:
        hashes = np.load(f"{p}.candidates.npy", mmap_mode="r")
        sample.append(np.asarray(hashes[::stride]))
        # ponytail: numpy probes, ~10^10 keys at full scale want a Rust kernel
        for i in range(0, len(hashes), BLOOM_CHUNK):
            pos = _positions(np.asarray(hashes[i : i + BLOOM_CHUNK]), len(bits) * 8).ravel()
            np.bitwise_or.at(
                bits, pos >> np.uint64(3), (1 << (pos & np.uint64(7))).astype(np.uint8)
            )
    t_max_hash = max(
        int(
            pl.scan_parquet(f"{p}.units.parquet")
            .select(pl.col("max_hash_g").max())
            .collect()
            .item()
            or 0
        )
        for p in prefixes
    )
    np.save(f"{out}.npy", bits)
    quantiles = (
        np.quantile(
            np.concatenate([np.zeros(0, dtype=np.uint64), *sample]),
            np.linspace(0, 1, QUANTILES + 1),
            method="inverted_cdf",
        )
        if n_keys
        else np.zeros(QUANTILES + 1, dtype=np.uint64)
    )
    meta = {
        "t_max_hash": t_max_hash,
        "keys": n_keys,
        "bits": len(bits) * 8,
        "quantiles": [int(q) for q in quantiles],
    }
    Path(f"{out}.json").write_text(json.dumps(meta) + "\n")


def presence(
    members_path: str | Path,
    prefix: str,
    bloom_prefix: str,
    bucket: int,
    n_ranges: int,
    params: IndexParams,
) -> None:
    """Stage 3: ``{prefix}.range{r}.parquet`` for each hash range r < ``n_ranges``, cut at
    quantiles of the candidate hashes so ranges hold about equal numbers of them.

    Rows: ``hash``, ``cluster_rep``, ``c``, ``candidate`` (hash <= the unit's own
    ``max_hash_g``) and ``bucket``, for every distinct (unit, k-mer) at ``t_max`` that
    passes the filter.
    """
    _check(params)
    meta = json.loads(Path(f"{bloom_prefix}.json").read_text())
    bits = np.load(f"{bloom_prefix}.npy")
    members, batches, _ = _prepare(members_path, params)
    units = _counting(members, pl.read_parquet(f"{prefix}.units.parquet"))
    rows = []
    for b in batches:
        kmers = _kmers(b, params, meta["quantiles"][-1])  # no candidate lies above it
        kmers = kmers.filter(bloom_contains(bits, kmers["hash"].to_numpy()))
        rows.append(kmers.group_by("unit", "hash").agg(c=pl.col("counts").sum().cast(pl.UInt32)))
    q = meta["quantiles"]
    cuts = np.array([q[r * (len(q) - 1) // n_ranges] for r in range(1, n_ranges)], np.uint64)
    table = (
        pl.concat(rows)
        .join(units.select("unit", "cluster_rep", "max_hash_g"), on="unit")
        .select(
            "hash",
            "cluster_rep",
            "c",
            candidate=pl.col("hash") <= pl.col("max_hash_g"),
            bucket=pl.lit(bucket, dtype=pl.UInt32),
        )
    )
    table = table.with_columns(
        range=pl.Series(
            np.searchsorted(cuts, table["hash"].to_numpy(), side="right"), dtype=pl.UInt32
        )
    )
    for r in range(n_ranges):
        table.filter(pl.col("range") == r).drop("range").write_parquet(f"{prefix}.range{r}.parquet")


def groups(paths: Sequence[str | Path], prefix: str) -> None:
    """Stage 4: one hash range's presence rows -> ``{prefix}.bucket{b}.parquet`` per bucket.

    Rows: candidate (``hash``, ``cluster_rep``) with ``c`` and ``n_groups``. Hashes with no
    candidate row anywhere are the filter's false positives and are dropped.
    """
    rows = pl.read_parquet([str(p) for p in paths])
    kept = rows.filter(pl.col("candidate").any().over("hash"))
    print(f"{prefix}: {rows.height - kept.height} of {rows.height} rows were false positives")
    scored = kept.with_columns(n_groups=pl.len().over("hash").cast(pl.UInt32)).filter("candidate")
    for (bucket,), part in scored.partition_by("bucket", as_dict=True).items():
        part.select("hash", "cluster_rep", "c", "n_groups").write_parquet(
            f"{prefix}.bucket{bucket}.parquet"
        )


def postings(
    members_path: str | Path,
    prefix: str,
    group_paths: Sequence[str | Path],
    params: IndexParams,
    pfam_path: str | Path | None = None,
) -> None:
    """Stage 5: ``{prefix}.postings.parquet``, ``.final.parquet`` (units with ``len_cv``),
    ``.pfam.parquet`` (with ``pfam_path``) and ``.final.json`` (stage-1 stats and this
    stage's; inputs stay unchanged, as Nextflow stages them as links)."""
    _check(params)
    members, batches, _ = _prepare(members_path, params)
    units = _counting(members, pl.read_parquet(f"{prefix}.units.parquet"))
    rows = (
        pl.read_parquet([str(p) for p in group_paths])
        if group_paths
        else pl.DataFrame(
            schema={
                "hash": pl.UInt64,
                "cluster_rep": units["cluster_rep"].dtype,
                "c": pl.UInt32,
                "n_groups": pl.UInt32,
            }
        )
    )
    rows = rows.join(units.select("cluster_rep", "unit"), on="cluster_rep").drop("cluster_rep")
    _, kept, stats = _select_postings(rows, units, params)
    units = units.join(
        _len_cv(batches, params, kept.select("unit", "hash"), members),
        on="unit",
        how="left",
        maintain_order="left",
    )
    rep = units.select("unit", "cluster_rep")
    _swap(kept, rep, "unit", "cluster_rep").write_parquet(f"{prefix}.postings.parquet")
    units.drop("unit").write_parquet(f"{prefix}.final.parquet")
    if pfam_path is not None:
        _swap(_unit_pfam(pfam_path, members), rep, "unit", "cluster_rep").write_parquet(
            f"{prefix}.pfam.parquet"
        )
    stats = json.loads(Path(f"{prefix}.stats.json").read_text()) | stats
    Path(f"{prefix}.final.json").write_text(json.dumps(stats) + "\n")


def pack(
    prefixes: Sequence[str],
    out_dir: str | Path,
    params: IndexParams,
    pfam: bool = False,
    postings_parquet: bool = False,
) -> dict[str, object]:
    """Stage 6: every bucket's units and postings -> the index in ``out_dir``."""
    _check(params)
    # ponytail: one process holds all postings; pack by hash range if this outgrows a node
    units = (
        pl.concat([pl.read_parquet(f"{p}.final.parquet") for p in prefixes], how="vertical_relaxed")
        .sort("cluster_rep")
        .with_columns(unit=pl.int_range(pl.len(), dtype=pl.UInt32))
    )
    rep = units.select("cluster_rep", "unit")

    def by_unit(name: str) -> pl.DataFrame:
        table = pl.concat([pl.read_parquet(f"{p}.{name}.parquet") for p in prefixes])
        return _swap(table, rep, "cluster_rep", "unit")

    stats: dict[str, object] = {}
    for p in prefixes:
        for key, value in json.loads(Path(f"{p}.final.json").read_text()).items():
            stats[key] = stats.get(key, 0) + value
    return finish_index(
        out_dir,
        params,
        units,
        by_unit("postings"),
        stats,
        unit_pfam=by_unit("pfam") if pfam else None,
        postings_parquet=postings_parquet,
    )
