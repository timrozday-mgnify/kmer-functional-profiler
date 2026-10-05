"""Partitioned index build (phase 6): the build of ``index.build_index`` split into jobs.

Members come in buckets that each hold whole clusters (``workflows/mgnify-subset``'s MERGE).
Everything a unit needs is local to its bucket except ``n_groups``, the number of units
holding a candidate k-mer, which feeds the score and the promiscuity cut. That is counted
exactly by a hash-partitioned reduce, with a Bloom filter of all candidate hashes to keep
the rows each bucket emits small. Units are keyed by ``cluster_rep`` until ``units``.

With a dense tier (``t_dense`` > 0) a unit's candidates are its k-mers with hash <= its
``max_hash_dense``, which contain its tier-2 candidates, so one reduce gives ``n_groups``
for both tiers; ``postings`` keeps the tier-2 ones as before, and the dense rows go through
``pack_range``, ``dedup`` and ``concat`` alongside tier 2 (files with a ``.dense`` infix).

1. ``candidates`` (per bucket): pass 1, the unit table and the bucket's candidate hashes.
2. ``bloom`` (once): a blocked Bloom filter (``_core.bloom_insert``) of every bucket's
   candidate hashes, and ``t_max``.
3. ``presence`` (per bucket): (hash, unit) of every k-mer at ``t_max`` that passes the
   filter, with its counting members and whether it is the unit's candidate, split into
   hash ranges.
4. ``groups`` (per hash range): ``n_groups`` per hash; keeps candidate rows only (which
   drops the filter's false positives) and splits them by bucket.
5. ``postings`` (per bucket): score, promiscuity cut and floor, as ``build_index``;
   the per-unit columns (``len_cv``, ``pin_hist``, ...) and Pfam labels.
6. ``units`` (once): the unit table, numbered in ``cluster_rep`` order, and the tier-2
   layout, with pack ranges cut at key boundaries.
7. ``pack_range`` (per pack range): that range's postings -> keys and value sets.
8. ``dedup`` (per set-hash range): every part's value sets in that range, deduplicated.
9. ``concat`` (once): the parts and deduplicated sets -> tier 2 and ``meta.json``.

The result equals ``build_index`` on the concatenated members.
"""

import json
import shutil
from collections.abc import Collection, Sequence
from pathlib import Path
from typing import Any, Final

import numpy as np
import polars as pl

from kmer_functional_profiler import _core
from kmer_functional_profiler.index import (
    PIN_BITS,
    PROMISCUOUS,
    IndexParams,
    PackedPart,
    PackedTable,
    SetSlice,
    _kmers,
    _len_cv,
    _max_hash_dense,
    _n_kmers,
    _prepare,
    _select_postings,
    _smallest,
    _unit_pfam,
    _unit_table,
    dedup_sets,
    dense_columns,
    dense_rows,
    key_shift,
    packed_layout,
    promiscuous_base,
    unit_columns,
    write_meta,
    write_unit_columns,
)

QUANTILES: Final = 1024  # candidate-hash quantiles recorded for balanced hash ranges
SAMPLE: Final = 1 << 20  # candidate hashes sampled for them
UNIT_COLUMNS: Final = (
    "cluster_rep",
    "n_members",
    "n_counting",
    "n_kmers",
    "t_g",
    "max_hash_g",
    "len_mean",
)
# The index's unit table, as build_index writes it (with ``unit`` first).
FINAL_COLUMNS: Final = (
    "cluster_rep",
    "n_members",
    "n_kmers",
    "t_g",
    "max_hash_g",
    "len_mean",
    "len_cv",
    "pin_hist",
    "pin_sum",
    "m_g",
)
# With a dense tier, as build_index writes them.
DENSE_FINAL_COLUMNS: Final = (
    *FINAL_COLUMNS[:7],
    "len_cv_dense",
    "max_hash_dense",
    *FINAL_COLUMNS[7:],
    "m_dense",
    "pin_hist_dense",
    "pin_sum_dense",
)


def _threshold(columns: Collection[str]) -> str:
    """The unit column bounding its candidates: the dense threshold if there is one."""
    return "max_hash_dense" if "max_hash_dense" in columns else "max_hash_g"


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
    members, batches, stats = _prepare(members_path, params)
    units = _unit_table(members, _n_kmers(batches, params), params)
    columns = list(UNIT_COLUMNS)
    if params.t_dense > 0:
        units = units.with_columns(max_hash_dense=_max_hash_dense(params))
        columns.append("max_hash_dense")
    bound = _threshold(columns)
    thresholds = units.select("unit", bound)
    t_max_hash = int(units[bound].max() or 0)  # type: ignore[arg-type]
    hashes = pl.concat(
        [
            _kmers(b, params, t_max_hash)
            .join(thresholds, on="unit")
            .filter(pl.col("hash") <= pl.col(bound))
            .select("hash")
            .unique()
            for b in batches
        ]
    ).unique()
    units.select(columns).write_parquet(f"{prefix}.units.parquet")
    np.save(f"{prefix}.candidates.npy", np.sort(hashes["hash"].to_numpy()))
    Path(f"{prefix}.stats.json").write_text(json.dumps(stats) + "\n")


def bloom(prefixes: Sequence[str], out: str | Path, bits_per_key: float = 10.0) -> None:
    """Stage 2: ``{out}.npy`` (filter bits) and ``{out}.json`` (``t_max_hash``, the largest
    candidate threshold, sizes and
    ``QUANTILES`` + 1 quantiles of the candidate hashes, where ``presence`` cuts ranges).

    Keys are counted per bucket, so hashes shared by buckets are counted more than once and
    the filter is at most that much larger than needed.
    """
    n_keys = sum(len(np.load(f"{p}.candidates.npy", mmap_mode="r")) for p in prefixes)
    block = _core.BLOOM_BLOCK_BYTES
    n_blocks = max(-(-int(n_keys * bits_per_key) // (8 * block)), 1)
    bits = np.zeros(n_blocks * block, dtype=np.uint8)
    stride, sample = max(n_keys // SAMPLE, 1), []
    for p in prefixes:
        hashes = np.load(f"{p}.candidates.npy", mmap_mode="r")
        sample.append(np.asarray(hashes[::stride]))
        _core.bloom_insert(bits, hashes)
    t_max_hash = max(
        int(
            pl.scan_parquet(f"{p}.units.parquet")
            .select(pl.col(_threshold(pl.read_parquet_schema(f"{p}.units.parquet"))).max())
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
    ``max_hash_g``, or ``max_hash_dense`` with a dense tier) and ``bucket``, for every
    distinct (unit, k-mer) at ``t_max_hash`` that passes the filter.
    """
    meta = json.loads(Path(f"{bloom_prefix}.json").read_text())
    bits = np.load(f"{bloom_prefix}.npy", mmap_mode="r")  # shared page cache across jobs
    members, batches, _ = _prepare(members_path, params)
    units = _counting(members, pl.read_parquet(f"{prefix}.units.parquet"))
    rows = []
    for b in batches:
        # Not the last quantile: that is the largest *sampled* candidate.
        kmers = _kmers(b, params, meta["t_max_hash"])
        kmers = kmers.filter(_core.bloom_contains(bits, kmers["hash"].to_numpy()))
        rows.append(kmers.group_by("unit", "hash").agg(c=pl.col("counts").sum().cast(pl.UInt32)))
    q = meta["quantiles"]
    cuts = np.array([q[r * (len(q) - 1) // n_ranges] for r in range(1, n_ranges)], np.uint64)
    bound = _threshold(units.columns)
    table = (
        pl.concat(rows)
        .join(units.select("unit", "cluster_rep", bound), on="unit")
        .select(
            "hash",
            "cluster_rep",
            "c",
            candidate=pl.col("hash") <= pl.col(bound),
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
    """Stage 5: ``{prefix}.postings.parquet`` (sorted by hash), ``.final.parquet`` (units
    with a posting and their index columns), ``.pfam.parquet`` (with ``pfam_path``),
    ``.dense.parquet`` (with a dense tier: those units' dense rows, sorted by hash),
    ``.promiscuous.npy`` (:func:`~kmer_functional_profiler.index.promiscuous_base`) and
    ``.final.json`` (stage-1 stats and this stage's; inputs stay unchanged, as Nextflow
    stages them as links)."""
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
    np.save(f"{prefix}.{PROMISCUOUS}", promiscuous_base(rows, params))
    units = units.join(
        _len_cv(batches, params, kept.select("unit", "hash"), members),
        on="unit",
        how="left",
        maintain_order="left",
    )
    rep = units.select("unit", "cluster_rep")
    _swap(kept, rep, "unit", "cluster_rep").sort("hash").write_parquet(f"{prefix}.postings.parquet")
    final, columns = unit_columns(units, kept), list(FINAL_COLUMNS)
    if params.t_dense > 0:
        dense = dense_rows(rows, units, params)
        final = dense_columns(
            final.join(
                _len_cv(batches, params, dense.select("unit", "hash"), members).rename(
                    {"len_cv": "len_cv_dense"}
                ),
                on="unit",
                how="left",
                maintain_order="left",
            ),
            dense,
        )
        posted = final.filter(pl.col("m_g") > 0).select("unit")
        _swap(dense.join(posted, on="unit", how="semi"), rep, "unit", "cluster_rep").sort(
            "hash"
        ).select("hash", "cluster_rep", "pin_q").write_parquet(f"{prefix}.dense.parquet")
        columns = list(DENSE_FINAL_COLUMNS)
    final.filter(pl.col("m_g") > 0).select(columns).write_parquet(f"{prefix}.final.parquet")
    if pfam_path is not None:
        _swap(_unit_pfam(pfam_path, members), rep, "unit", "cluster_rep").write_parquet(
            f"{prefix}.pfam.parquet"
        )
    stats = json.loads(Path(f"{prefix}.stats.json").read_text()) | stats
    Path(f"{prefix}.final.json").write_text(json.dumps(stats) + "\n")


def units(
    prefixes: Sequence[str],
    bloom_prefix: str,
    n_ranges: int,
    params: IndexParams,
    out_dir: str | Path,
    pfam: bool = False,
) -> None:
    """Stage 6: ``units.parquet``, ``unit_pfam.parquet`` (with ``pfam``),
    ``promiscuous.npy`` (the buckets' merged) and ``pack.json``
    (tier-2 layout, pack range bounds, stats, and with a dense tier its layout and bounds
    under ``dense``) in ``out_dir``.

    Pack ranges are cut at the candidate-hash quantiles, rounded down to key boundaries so
    that hashes sharing a key never fall in two ranges.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    finals = [f"{p}.final.parquet" for p in prefixes]
    dense = params.t_dense > 0
    # Each bucket sized its p_in histograms by its own largest count; use the overall one.
    hists = {"pin_hist": "m_g", **({"pin_hist_dense": "m_dense"} if dense else {})}
    casts = []
    for name, count in hists.items():
        max_m = pl.scan_parquet(finals).select(pl.col(count).max()).collect().item() or 0
        casts.append(
            pl.col(name).cast(pl.Array(pl.Series(np.zeros(0, _smallest(max_m))).dtype, 2**PIN_BITS))
        )
    table = (
        pl.concat([pl.read_parquet(f).with_columns(casts) for f in finals])
        .sort("cluster_rep")
        .select(pl.int_range(pl.len(), dtype=pl.UInt32).alias("unit"), pl.all())
    )
    table.write_parquet(out / "units.parquet")
    rep = table.select("cluster_rep", "unit")
    if pfam:
        pfams = pl.concat([pl.read_parquet(f"{p}.pfam.parquet") for p in prefixes])
        _swap(pfams, rep, "cluster_rep", "unit").sort("unit", "pfam_accession").write_parquet(
            out / "unit_pfam.parquet"
        )
    promiscuous = np.unique(np.concatenate([np.load(f"{p}.{PROMISCUOUS}") for p in prefixes]))
    np.save(out / PROMISCUOUS, promiscuous)
    stats: dict[str, object] = {"promiscuous_base": len(promiscuous)}
    for p in prefixes:
        for key, value in json.loads(Path(f"{p}.final.json").read_text()).items():
            stats[key] = stats.get(key, 0) + value
    n_postings = int(table["m_g"].sum())
    stats |= {
        "n_units": table.height,
        "t_max": float(table["t_g"].max() or 0.0),  # type: ignore[arg-type]
        "postings": n_postings,
    }
    q = json.loads(Path(f"{bloom_prefix}.json").read_text())["quantiles"]

    def pack(max_hash: str, n: int) -> dict[str, object]:
        layout = packed_layout(int(table[max_hash].max() or 0), params.fp_bits, n)  # type: ignore[arg-type]
        shift = key_shift(layout)
        cuts = [q[r * (len(q) - 1) // n_ranges] >> shift << shift for r in range(1, n_ranges)]
        return {"layout": layout, "bounds": [0, *cuts, 2**64]}

    meta = pack("max_hash_g", n_postings) | {"stats": stats}
    if dense:
        meta["dense"] = pack("max_hash_dense", int(table["m_dense"].sum()))
    (out / "pack.json").write_text(json.dumps(meta) + "\n")


def pack_range(prefixes: Sequence[str], units_dir: str | Path, part: int, out_prefix: str) -> None:
    """Stage 7: postings with hash in pack range ``part`` -> ``PackedPart`` files at
    ``out_prefix`` and ``{out_prefix}.json`` (postings and distinct hashes); with a dense
    tier, its rows in its own range ``part`` likewise at ``{out_prefix}.dense``."""
    meta = json.loads((Path(units_dir) / "pack.json").read_text())
    _pack([f"{p}.postings.parquet" for p in prefixes], units_dir, meta, part, out_prefix)
    if "dense" in meta:
        files = [f"{p}.dense.parquet" for p in prefixes]
        _pack(files, units_dir, meta["dense"], part, f"{out_prefix}.dense")


def _pack(
    files: list[str], units_dir: str | Path, meta: dict[str, Any], part: int, out_prefix: str
) -> None:
    """Rows of ``files`` in pack range ``part`` of ``meta`` (layout and bounds), packed."""
    lo, hi = meta["bounds"][part], meta["bounds"][part + 1]
    in_range = pl.col("hash") >= lo
    if hi < 2**64:
        in_range &= pl.col("hash") < hi
    rows = (
        pl.scan_parquet(files)
        .filter(in_range)
        .select("hash", "cluster_rep", "pin_q")
        .join(
            pl.scan_parquet(Path(units_dir) / "units.parquet").select("cluster_rep", "unit"),
            on="cluster_rep",
        )
        .collect()
    )
    values = (rows["unit"].cast(pl.UInt64).to_numpy() << np.uint64(PIN_BITS)) | rows[
        "pin_q"
    ].to_numpy()
    PackedPart.build(rows["hash"].to_numpy(), values, meta["layout"]).save(out_prefix)
    counts = {"postings": rows.height, "distinct_hashes": rows["hash"].n_unique()}
    Path(f"{out_prefix}.json").write_text(json.dumps(counts) + "\n")


def dedup(part_prefixes: Sequence[str], index: int, n_ranges: int, out_prefix: str) -> None:
    """Stage 8: the parts' value sets whose first hash word lies in set-hash range
    ``index`` of ``n_ranges`` (equal widths; the hashes are uniform) -> ``SetSlice`` files
    at ``out_prefix``; the dense parts' sets, if any, at ``{out_prefix}.dense``."""
    lo = index * 2**64 // n_ranges
    hi = None if index == n_ranges - 1 else (index + 1) * 2**64 // n_ranges
    dedup_sets([PackedPart.load(p) for p in part_prefixes], lo, hi).save(out_prefix)
    if part_prefixes and Path(f"{part_prefixes[0]}.dense.json").exists():
        dense = [PackedPart.load(f"{p}.dense") for p in part_prefixes]
        dedup_sets(dense, lo, hi).save(f"{out_prefix}.dense")


def concat(
    part_prefixes: Sequence[str],
    set_prefixes: Sequence[str],
    units_dir: str | Path,
    out_dir: str | Path,
    params: IndexParams,
) -> dict[str, object]:
    """Stage 9: pack range parts and set slices (each in range order) -> tier 2, the unit
    tables and ``meta.json`` in ``out_dir``; returns the stats."""
    units_dir, out = Path(units_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    meta = json.loads((units_dir / "pack.json").read_text())
    counts = [json.loads(Path(f"{p}.json").read_text()) for p in part_prefixes]
    if sum(c["postings"] for c in counts) != meta["stats"]["postings"]:
        raise ValueError("pack ranges do not hold every posting")
    tables = {
        "tier2": PackedTable.concat(
            meta["layout"],
            [PackedPart.load(p) for p in part_prefixes],
            [SetSlice.load(p) for p in set_prefixes],
        )
    }
    stats = meta["stats"] | {"distinct_hashes": sum(c["distinct_hashes"] for c in counts)}
    if "dense" in meta:
        dense = [json.loads(Path(f"{p}.dense.json").read_text())["postings"] for p in part_prefixes]
        tables["dense"] = PackedTable.concat(
            meta["dense"]["layout"],
            [PackedPart.load(f"{p}.dense") for p in part_prefixes],
            [SetSlice.load(f"{p}.dense") for p in set_prefixes],
        )
        stats |= {"dense_postings": sum(dense), "dense_bytes": tables["dense"].nbytes()}
    for name in ("units.parquet", "unit_pfam.parquet", PROMISCUOUS):
        if (units_dir / name).exists() and units_dir.resolve() != out.resolve():
            shutil.copyfile(units_dir / name, out / name)
    write_unit_columns(out)
    return write_meta(out, params, tables, stats)
