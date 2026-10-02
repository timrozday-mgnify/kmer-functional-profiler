"""Mask sidecar (Additional references, step 2): an index's k-mers found in a host genome.

``build_mask`` six-frame translates a genome (FASTA, plain or gzip) in overlapping windows
through the read kernel, keeps the hashes the index holds (tier 2 or dense), and writes

- ``mask.npy``: those hashes, sorted; the query drops sampled hashes in it before lookup;
- ``mask_units.parquet`` (and ``mask_units_dense.parquet``): per unit of the index with
  masked postings in tier 2 (dense), their count, per-``p_in``-level counts and ``p_in`` sum,
  which the query subtracts from ``m_g``, ``pin_hist``, ``pin_sum`` (``m_dense``, ...) so
  expected counts match the k-mers that can still be hit;
- ``mask.json``: the genome's path and SHA-256, the index's ``meta.json`` SHA-256, counts.

Found hashes come from fingerprint lookups, so a host hash that only shares a fingerprint
with an index k-mer is masked too: harmless for the query (masking compares full hashes),
but it decrements its set's units by a k-mer they keep (~6e-5 per host lookup at 16
fingerprint bits). ``p_in`` decrements use the quantised level, as the tables store it.
"""

import gzip
import hashlib
import json
from collections.abc import Iterator
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl

from kmer_functional_profiler import _core
from kmer_functional_profiler.index import PIN_BITS, Index
from kmer_functional_profiler.query import unit_hits

WINDOW: Final = 10_000  # nt per window
WINDOWS_PER_BATCH: Final = 2_000  # windows hashed per kernel call (~20 Mbp)
HOST_LIKE: Final = 0.5  # masked fraction of a unit's tier-2 k-mers above which it is host_like
LEVELS: Final = 2**PIN_BITS


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def windows(path: Path, size: int, overlap: int) -> Iterator[bytes]:
    """Each FASTA record cut into windows of ``size`` nt overlapping by ``overlap``, so
    every stretch of ``overlap + 1`` nt lies whole in some window."""
    opener = gzip.open if path.suffix == ".gz" else open
    buf = bytearray()
    with opener(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if len(buf) > overlap or (buf and not overlap):
                    yield bytes(buf)
                buf.clear()
                continue
            buf += line.strip()
            while len(buf) >= size:
                yield bytes(buf[:size])
                del buf[: size - overlap]
    if len(buf) > overlap:
        yield bytes(buf)


def _batches(path: Path, k: int) -> Iterator[list[bytes]]:
    batch: list[bytes] = []
    for w in windows(path, WINDOW, 3 * k - 1):
        batch.append(w)
        if len(batch) == WINDOWS_PER_BATCH:
            yield batch
            batch = []
    if batch:
        yield batch


def _decrements(index: Index, tier: str, hashes: np.ndarray) -> pl.DataFrame:
    """Per unit: masked postings of ``tier``, per p_in level and their p_in sum."""
    table = getattr(index, tier)
    column = "max_hash_g" if tier == "tier2" else "max_hash_dense"
    hashes = hashes[hashes <= np.uint64(table.max_hash)]
    hits = unit_hits(table, index.units[column], hashes, np.zeros(len(hashes), np.uint64))
    unit = hits["unit"].to_numpy().astype(np.int64)
    level = hits["pin_q"].to_numpy().astype(np.int64)
    units, at = np.unique(unit, return_inverse=True)
    hist = np.zeros((len(units), LEVELS), dtype=np.uint32)
    np.add.at(hist, (at, level), 1)
    return pl.DataFrame(
        {
            "unit": units.astype(np.uint32),
            "masked": hist.sum(axis=1).astype(np.uint32),
            "pin_hist": hist,
            "pin_sum": hist @ (np.arange(LEVELS) / (LEVELS - 1)),
        }
    )


def build_mask(genome: str | Path, index_dir: str | Path, out: str | Path) -> dict[str, object]:
    """Write the mask sidecar of ``genome`` for the index in ``index_dir`` to ``out``."""
    genome, index_dir, out = Path(genome), Path(index_dir), Path(out)
    index = Index.load(index_dir)
    params = index.meta["params"]
    if index.meta.get("hash") == "sourmash":
        raise ValueError("masks need the index's own hash (not a sourmash import)")
    tiers = ["tier2"] + (["dense"] if index.dense is not None else [])
    max_hash = max(int(getattr(index, t).max_hash) for t in tiers)
    found, sampled = [], 0
    for batch in _batches(genome, params["k"]):
        raw = _core.hash_dna(
            batch, params["k"], alphabet=params["alphabet"], frames="all", max_hash=max_hash
        )["hash"]
        sampled += len(raw)
        hashes = np.unique(raw)
        held = np.zeros(len(hashes), dtype=bool)
        for t in tiers:
            table = getattr(index, t)
            ok = hashes <= np.uint64(table.max_hash)
            held[ok] |= table.lookup(hashes[ok]) >= 0
        found.append(hashes[held])
    mask = np.unique(np.concatenate([np.empty(0, np.uint64), *found]))
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "mask.npy", mask)
    units = _decrements(index, "tier2", mask)
    units.write_parquet(out / "mask_units.parquet")
    if index.dense is not None:
        _decrements(index, "dense", mask).write_parquet(out / "mask_units_dense.parquet")
    stats = {
        "genome": str(genome),
        "genome_sha256": _sha256(genome),
        "index_meta_sha256": _sha256(index_dir / "meta.json"),
        "sampled_kmers": sampled,
        "masked_hashes": len(mask),
        "masked_postings": int(units["masked"].sum()),
        "postings": index.meta["stats"].get("postings"),
        "units_masked": units.height,
    }
    (out / "mask.json").write_text(json.dumps(stats, indent=2) + "\n")
    return stats


class Mask:
    """A loaded mask sidecar: drops masked hashes and adjusts the unit rows of its index."""

    def __init__(self, directory: str | Path) -> None:
        directory = Path(directory)
        self.hashes: np.ndarray = np.load(directory / "mask.npy")
        self.meta = json.loads((directory / "mask.json").read_text())
        self.units = {
            suffix: pl.read_parquet(path)
            for suffix, path in (
                ("", directory / "mask_units.parquet"),
                ("_dense", directory / "mask_units_dense.parquet"),
            )
            if path.exists()
        }

    def keep(self, hashes: np.ndarray) -> np.ndarray:
        """Whether each hash is not masked."""
        if not len(self.hashes):
            return np.ones(len(hashes), dtype=bool)
        at = np.minimum(np.searchsorted(self.hashes, hashes), len(self.hashes) - 1)
        out: np.ndarray = self.hashes[at] != hashes
        return out

    def adjust(self, rows: pl.DataFrame) -> pl.DataFrame:
        """``rows`` of the unit table (``unit``: the index's ids, sorted) less their masked
        postings, with ``masked_fraction`` (of tier-2 k-mers) and ``host_like`` (above
        ``HOST_LIKE``). Counts stay >= 1: a unit masked whole is hit only by fingerprint
        false hits."""
        unit = rows["unit"].to_numpy().astype(np.int64)
        out: dict[str, np.ndarray] = {}
        for suffix, dec in self.units.items():
            m = "m_g" if suffix == "" else "m_dense"
            if m not in rows.columns:
                continue
            ids = dec["unit"].to_numpy().astype(np.int64)
            at = np.minimum(np.searchsorted(ids, unit), max(len(ids) - 1, 0))
            hit = ids[at] == unit if len(ids) else np.zeros(len(unit), dtype=bool)
            masked = np.zeros(len(unit))
            masked[hit] = dec["masked"].to_numpy()[at[hit]]
            if suffix == "":
                out["masked_fraction"] = masked / rows["m_g"].to_numpy()
            out[m] = np.maximum(rows[m].to_numpy() - masked, 1)
            pin_sum = rows[f"pin_sum{suffix}"].to_numpy().astype(np.float64)
            pin_sum[hit] -= dec["pin_sum"].to_numpy()[at[hit]]
            out[f"pin_sum{suffix}"] = np.maximum(pin_sum, 1e-3)
            hist = rows[f"pin_hist{suffix}"].to_numpy().astype(np.int64)
            hist[hit] -= dec["pin_hist"].to_numpy()[at[hit]].astype(np.int64)
            out[f"pin_hist{suffix}"] = np.maximum(hist, 0)
        fraction = out.pop("masked_fraction", np.zeros(len(unit)))
        return rows.with_columns(
            *(pl.Series(c, v).cast(rows.schema[c]) for c, v in out.items()),
            masked_fraction=pl.Series(fraction),
            host_like=pl.Series(fraction > HOST_LIKE),
        )
