"""Query (phase 3 prototype): naive hit counts and containment per unit.

Reads are streamed through the Rust kernel at the index's ``t_max``; every sampled hash
is looked up in tier 2 and counts for each unit in its set whose own threshold it passes
(``hash <= max_hash_g``), which also discards most fingerprint false hits. No query
sketch is stored: only per-(unit, hash) counts are kept. Indexes imported from sourmash
signatures (``meta["hash"] == "sourmash"``) hash reads with sourmash instead; ``frames``
and ``genetic_code`` then do not apply.

Shared k-mers count for every unit that holds them in ``kmers_hit``; ``kmers_unique`` is
what a gather-style greedy assignment leaves each unit (the phase-4 detection baseline).
"""

import heapq
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import polars as pl

from kmer_functional_profiler import _core
from kmer_functional_profiler.compat import sourmash_hits
from kmer_functional_profiler.index import PIN_BITS, Index, IndexParams


def unit_hits(index: Index, hashes: np.ndarray, reads: np.ndarray) -> pl.DataFrame:
    """Expand sampled query hashes to (unit, hash, read) hits through tier 2."""
    tier2 = index.tier2
    set_ids = tier2.lookup(hashes)
    found = set_ids >= 0
    hashes, reads, set_ids = hashes[found], reads[found], set_ids[found]
    starts = tier2.set_offsets[set_ids].astype(np.int64)
    lengths = tier2.set_offsets[set_ids + 1].astype(np.int64) - starts
    # Positions of every set member, hit by hit.
    within = np.arange(lengths.sum()) - np.repeat(np.cumsum(lengths) - lengths, lengths)
    values = tier2.set_values[np.repeat(starts, lengths) + within]
    units = (values >> np.uint64(PIN_BITS)).astype(np.uint32)
    hashes, reads = np.repeat(hashes, lengths), np.repeat(reads, lengths)
    keep = hashes <= index.units["max_hash_g"].to_numpy()[units]
    return pl.DataFrame({"unit": units[keep], "hash": hashes[keep], "read": reads[keep]})


def gather(kmers: pl.DataFrame, t_g: np.ndarray) -> pl.DataFrame:
    """Assign each hit k-mer to one unit, greedily, as ``sourmash gather`` does.

    Repeatedly take the unit with the most unassigned hit k-mers, scaled by ``1 / t_g`` so
    units sampled at different rates compare as estimated k-mer counts, and give it those
    k-mers. Units left with none are explained away and get no row. Returns ``unit``,
    ``kmers_unique`` (k-mers assigned) and ``gather_rank`` (0 = taken first).
    """
    remaining = {u: set(h) for u, h in kmers.group_by("unit").agg("hash").iter_rows()}
    # Scores only fall as k-mers are taken, so a lazy heap needs only stale-top rechecks.
    heap = [(-len(h) / t_g[u], u) for u, h in remaining.items()]
    heapq.heapify(heap)
    taken: set[int] = set()
    rows: list[tuple[int, int, int]] = []
    while heap:
        _, unit = heapq.heappop(heap)
        mine = remaining[unit] = remaining[unit] - taken  # iterates the unit's set, not taken
        if not mine:
            continue
        score = len(mine) / t_g[unit]
        if heap and score < -heap[0][0]:
            heapq.heappush(heap, (-score, unit))
            continue
        taken |= mine
        rows.append((unit, len(mine), len(rows)))
    return pl.DataFrame(
        rows,
        schema={"unit": pl.UInt32, "kmers_unique": pl.UInt32, "gather_rank": pl.UInt32},
        orient="row",
    )


def profile(
    index: Index,
    r1: str | Path,
    r2: str | Path | None = None,
    *,
    genetic_code: int = 11,
    frames: str = "stopfree",
    batch_reads: int = 100_000,
) -> pl.DataFrame:
    """Per-unit hits, distinct k-mers hit, reads hit, containment and mean coverage.

    ``containment`` is the fraction of the unit's kept k-mers seen at least once;
    ``coverage`` is hits per kept k-mer; ``kmers_unique`` and ``gather_rank`` come from
    :func:`gather` (0 and null for units explained away). Units without hits are omitted.
    """
    params = IndexParams(**index.meta["params"])
    batches: Iterable[dict[str, np.ndarray]]
    if index.meta.get("hash") == "sourmash":
        batches = sourmash_hits(r1, r2, params.k, index.tier2.max_hash, batch_reads)
    else:
        batches = _core.FastxHits(
            r1,
            r2,
            k=params.k,
            alphabet=params.alphabet,
            genetic_code=genetic_code,
            frames=frames,
            max_hash=index.tier2.max_hash,
            batch_reads=batch_reads,
        )

    def counts(hashes: np.ndarray, reads: np.ndarray) -> pl.DataFrame:
        hits = unit_hits(index, hashes, reads)
        return hits.group_by("unit", "hash").agg(hits=pl.len(), reads=pl.col("read").unique())

    empty = np.empty(0, dtype=np.uint64)
    per_kmer = pl.concat([counts(empty, empty), *(counts(b["hash"], b["read"]) for b in batches)])
    assigned = gather(per_kmer.select("unit", "hash").unique(), index.units["t_g"].to_numpy())
    return (
        per_kmer.group_by("unit")
        .agg(
            hits=pl.col("hits").sum().cast(pl.UInt64),
            kmers_hit=pl.col("hash").n_unique().cast(pl.UInt32),
            reads=pl.col("reads").explode(empty_as_null=False).n_unique().cast(pl.UInt64),
        )
        .join(index.units, on="unit")
        .join(assigned, on="unit", how="left")
        .with_columns(
            containment=pl.col("kmers_hit") / pl.col("m_g"),
            coverage=pl.col("hits") / pl.col("m_g"),
            kmers_unique=pl.col("kmers_unique").fill_null(0),
        )
        .sort("unit")
    )
