"""Partitioned build: equal to the single-process build; Bloom filter has no false negatives."""

import itertools
import json
import random
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from kmer_functional_profiler import _core, partition
from kmer_functional_profiler.index import (
    IndexParams,
    PackedPart,
    PackedTable,
    build_index,
    key_shift,
    packed_layout,
)

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
BUCKETS, RANGES = 3, 4


def test_bloom_filter() -> None:
    rng = np.random.default_rng(1)
    keys = rng.integers(0, 2**62, 20_000, dtype=np.uint64)
    bits = np.zeros(20_000 * 10 // 8 // 64 * 64, dtype=np.uint8)
    _core.bloom_insert(bits, keys)
    assert _core.bloom_contains(bits, keys).all()
    others = rng.integers(0, 2**62, 100_000, dtype=np.uint64)
    assert _core.bloom_contains(bits, others).mean() < 0.02  # ~1% expected
    with pytest.raises(ValueError, match="multiple of 64"):
        _core.bloom_insert(np.zeros(65, dtype=np.uint8), keys)


def test_packed_table_from_parts_equals_whole() -> None:
    rng = np.random.default_rng(2)
    max_hash = 2**40 - 1
    hashes = rng.integers(0, max_hash, 3000, dtype=np.uint64)
    values = rng.integers(0, 50, 3000, dtype=np.uint64)  # repeated values -> shared sets
    # 4 fingerprint bits: many keys collide and pool their values.
    whole = PackedTable.build(hashes, values, max_hash, 4)
    layout = packed_layout(max_hash, 4, len(hashes))
    shift = key_shift(layout)
    cuts = [0, *(int(c) >> shift << shift for c in np.quantile(hashes, [0.3, 0.6])), 2**64]
    parts = []
    for lo, hi in itertools.pairwise(cuts):
        in_range = (hashes >= np.uint64(lo)) & (hashes < hi) if hi < 2**64 else hashes >= lo
        parts.append(PackedPart.build(hashes[in_range], values[in_range], layout))
    joined = PackedTable.concat(layout, parts)
    for field in ("offsets", "fingerprints", "set_ids", "set_offsets", "set_values"):
        a, b = getattr(whole, field), getattr(joined, field)
        assert a.dtype == b.dtype and np.array_equal(a, b), field
    assert len(whole.set_offsets) - 1 < len(whole.fingerprints)  # keys share sets


def test_partitioned_build_equals_single(tmp_path: Path) -> None:
    rng = random.Random(3)
    rows, pfam = [], []
    for rep in range(1, 300):
        base = "".join(rng.choices(AMINO_ACIDS, k=rng.randint(20, 300)))
        for m in range(rng.choice([1, 1, 2, 3, 8])):
            seq = "".join(rng.choice(AMINO_ACIDS) if rng.random() < 0.05 else a for a in base)
            rows.append((rep * 100 + m, rep * 100, rng.random() < 0.7, seq))
            pfam += [(rep * 100 + m, f"PF{rep % 7:05d}")] * (rng.random() < 0.5)
    members = pl.DataFrame(
        rows, schema=["protein_id", "cluster_rep", "full_length", "sequence"], orient="row"
    )
    members.write_parquet(tmp_path / "members.parquet")
    pl.DataFrame(pfam, schema=["protein_id", "pfam_accession"], orient="row").write_parquet(
        tmp_path / "pfam.parquet"
    )
    # k = 4 shares many k-mers across buckets, so n_groups and the promiscuity cut matter.
    params = IndexParams(k=4, t_base=0.05, n_min=4, t_cap=0.5, max_groups=3)
    single = build_index(
        tmp_path / "members.parquet", tmp_path / "single", params, tmp_path / "pfam.parquet", True
    )
    assert single["promiscuous_dropped"] > 0  # type: ignore[operator]

    d = tmp_path / "parts"
    d.mkdir()
    prefixes = [str(d / f"b{b}") for b in range(BUCKETS)]
    for b, prefix in enumerate(prefixes):
        members.filter(pl.col("cluster_rep") % BUCKETS == b).write_parquet(
            f"{prefix}.members.parquet"
        )
        partition.candidates(f"{prefix}.members.parquet", prefix, params)
    # A small filter, so many rows are false positives that groups must drop.
    partition.bloom(prefixes, d / "bloom", bits_per_key=2)
    for b, prefix in enumerate(prefixes):
        partition.presence(f"{prefix}.members.parquet", prefix, str(d / "bloom"), b, RANGES, params)
    for r in range(RANGES):
        partition.groups([f"{p}.range{r}.parquet" for p in prefixes], str(d / f"range{r}"))
    for b, prefix in enumerate(prefixes):
        found = [d / f"range{r}.bucket{b}.parquet" for r in range(RANGES)]
        partition.postings(
            f"{prefix}.members.parquet",
            prefix,
            [p for p in found if p.exists()],
            params,
            tmp_path / "pfam.parquet",
        )
    partition.units(prefixes, str(d / "bloom"), RANGES, params, d / "units", pfam=True)
    parts = [str(d / f"part{r}") for r in range(RANGES)]
    for r, part in enumerate(parts):
        partition.pack_range(prefixes, d / "units", r, part)
    stats = partition.concat(parts, d / "units", tmp_path / "parted", params)

    assert stats.keys() == single.keys()
    for key, value in single.items():
        assert stats[key] == pytest.approx(value), key
    for name in ("units", "unit_pfam"):
        a, b = (pl.read_parquet(tmp_path / x / f"{name}.parquet") for x in ("single", "parted"))
        assert a.equals(b), name
    for npy in (tmp_path / "single").glob("*.npy"):
        assert np.array_equal(np.load(npy), np.load(tmp_path / "parted" / npy.name)), npy.name
    meta = [json.loads((tmp_path / x / "meta.json").read_text()) for x in ("single", "parted")]
    assert meta[0]["tier2"] == meta[1]["tier2"]
