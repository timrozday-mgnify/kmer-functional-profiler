"""sourmash compatibility: protein FracMinHash signatures as an index, sourmash read hashing.

``import_signatures`` turns a signature collection (e.g. fmh-funprofiler's KO sketches)
into an index with one unit per signature, ``p_in = 1`` and no promiscuity cut, marked
``"hash": "sourmash"`` in ``meta.json``. Querying such an index hashes reads the way
``sourmash sketch translate`` does, so per-unit ``kmers_hit`` equals the overlap
``sourmash prefetch`` reports and ``containment`` its ``f_match``.
"""

from collections.abc import Iterator
from pathlib import Path

import numpy as np
import polars as pl
import screed
import sourmash
from numpy.typing import NDArray

from kmer_functional_profiler.index import PIN_BITS, IndexParams, write_index


def import_signatures(sig_path: str | Path, out_dir: str | Path, ksize: int) -> dict[str, object]:
    """Write an index of the protein signatures with ``ksize`` in ``sig_path``."""
    collection = sourmash.load_file_as_index(str(sig_path)).select(ksize=ksize, moltype="protein")
    names, sizes, hashes = [], [], []
    scaled, max_hash = None, None
    for sig in collection.signatures():
        mh = sig.minhash
        scaled, max_hash = scaled or mh.scaled, max_hash or mh._max_hash
        if (mh.scaled, mh._max_hash) != (scaled, max_hash):
            raise ValueError(f"{sig.name}: scaled {mh.scaled} differs from {scaled}")
        names.append(sig.name)
        sizes.append(len(mh))
        hashes.append(np.fromiter(mh.hashes, dtype=np.uint64, count=len(mh)))
    if scaled is None or max_hash is None:
        raise ValueError(f"no protein signatures with ksize {ksize} in {sig_path}")

    units = pl.DataFrame(
        {"unit": np.arange(len(names), dtype=np.uint32), "name": names},
    ).with_columns(t_g=pl.lit(1 / scaled), max_hash_g=pl.lit(max_hash, dtype=pl.UInt64))
    scored = (
        pl.DataFrame(
            {
                "hash": np.concatenate(hashes),
                "unit": np.repeat(units["unit"].to_numpy(), sizes),
            }
        )
        .with_columns(n_groups=pl.len().over("hash").cast(pl.UInt32))
        .with_columns(
            p_in=pl.lit(1.0),
            pin_q=pl.lit(2**PIN_BITS - 1, dtype=pl.UInt64),
            score=-pl.col("n_groups").log(2),
        )
    )
    params = IndexParams(k=ksize, t_base=1 / scaled, n_min=0, max_groups=2**32)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    stats = {"source": str(sig_path), "scaled": scaled}
    return write_index(out, params, units, scored, stats, hash_scheme="sourmash")


def sourmash_hits(
    r1: str | Path,
    r2: str | Path | None,
    ksize: int,
    max_hash: int,
    batch_reads: int = 100_000,
) -> Iterator[dict[str, NDArray[np.uint64]]]:
    """Batches of ``hash``/``read`` columns hashed like ``sourmash sketch translate``.

    Both mates of a pair share a read index. Hashes above ``max_hash`` are dropped.
    """
    mh = sourmash.MinHash(n=0, ksize=ksize, is_protein=True, scaled=1)
    hashes: list[NDArray[np.uint64]] = []
    reads: list[NDArray[np.uint64]] = []
    for read, pair in enumerate(_pairs(r1, r2)):
        h = np.array(
            [x for rec in pair for x in mh.seq_to_hashes(rec.sequence, force=True)],
            dtype=np.uint64,
        )
        h = h[h <= np.uint64(max_hash)]
        hashes.append(h)
        reads.append(np.full(len(h), read, dtype=np.uint64))
        if (read + 1) % batch_reads == 0:
            yield {"hash": np.concatenate(hashes), "read": np.concatenate(reads)}
            hashes, reads = [], []
    if hashes:
        yield {"hash": np.concatenate(hashes), "read": np.concatenate(reads)}


def _pairs(r1: str | Path, r2: str | Path | None) -> Iterator[tuple[screed.Record, ...]]:
    with screed.open(str(r1)) as f1:
        if r2 is None:
            yield from ((rec,) for rec in f1)
            return
        with screed.open(str(r2)) as f2:
            yield from zip(f1, f2, strict=True)
