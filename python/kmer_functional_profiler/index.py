"""Index build (phase 2 prototype): MGnify90 clusters as units, sampled member k-mers.

Input is a members table (``protein_id``, ``cluster_rep``, ``full_length``, ``sequence``),
e.g. from ``workflows/mgnify-subset``. Translated Illumina adapters in member sequences are
masked first (``mask_adapters``). The build then hashes every member three times with the
Rust kernel, keeping only small tables in memory:

1. all k-mers -> distinct k-mers per unit, ``n_kmers``;
2. at ``t_max`` -> candidate hashes, i.e. (unit, k-mer) with hash <= the unit's threshold
   ``t_g = max(t_base, min(t_cap, oversample * n_min / n_kmers))`` (``t_base`` for
   singletons); ``t_cap`` bounds ``t_max`` and so the query's sampling rate;
3. at ``t_max`` -> presence of candidate hashes in every unit, giving ``p_in`` (fraction of
   full-length members containing the k-mer), ``n_groups`` and the score
   ``log2(p_in / n_groups)``.

Postings are candidates minus promiscuous k-mers (``n_groups > max_groups``); units whose
``t_g`` was raised above ``t_base`` keep only their ``n_min`` best-scoring candidates. The
query stays consistent because any subset of k-mers with hash <= ``t_g`` may be kept. Tier 2
maps each posting hash to its set of (unit, quantised ``p_in``), stored as bucketed
fingerprints (``PackedTable``) in ``.npy`` files. Only units with a posting can be hit, so
only they get a row in ``units.parquet``, renumbered from 0 in ``cluster_rep`` order;
``postings.parquet`` (every posting with its score) is written only for inspection.

With ``t_dense`` > 0 a fourth pass keeps every unit's k-mers with hash <= max(``t_dense``,
``t_g``) (``max_hash_dense``, so each unit's dense set contains its tier-2 set), minus
promiscuous ones, in a ``dense`` table of the same layout; ``m_dense`` counts them per unit.
``pin_hist`` (``pin_hist_dense``) counts each unit's kept (dense) k-mers per ``p_in`` level,
and ``pin_sum`` (``pin_sum_dense``) sums their ``p_in``: how many of them an average member
holds, which turns present k-mers into member-equivalents (copies). ``len_cv``
(``len_cv_dense``) is the coefficient of variation of kept k-mers over the counting
members: how far one copy's k-mers can stray from ``pin_sum``.
The query probes it only for the units the sparse tier detects, to fit abundances on more k-mers.
"""

import json
from collections.abc import Iterator, Sequence
from dataclasses import asdict, dataclass
from functools import cached_property
from pathlib import Path
from typing import Any, Final, NamedTuple, Self

import duckdb
import numpy as np
import polars as pl
from numpy.typing import NDArray

from kmer_functional_profiler import _core

U64_MAX: Final = 2**64 - 1
PIN_BITS: Final = 4  # tier-2 values are unit << PIN_BITS | quantised p_in
PACKED_FIELDS: Final = ("offsets", "fingerprints", "set_ids", "set_offsets", "set_values")
KEYS_PER_BUCKET_BITS: Final = 2  # ~4 keys per offsets bucket
SET_SEEDS: Final = (0x243F6A8885A308D3, 0x13198A2E03707344)  # two 64-bit set-content hashes
# Constant parts of Illumina TruSeq and Nextera adapters and the P5/P7 flow-cell ends. Their
# translations recur in MGnify proteins predicted from reads with adapter read-through.
ADAPTERS: Final = (
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
    "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
    "ATCTCGTATGCCGTCTTCTGCTTG",
    "GTGTAGATCTCGGTGGTCGCCGTATCATT",
)
MASK_WIDTH: Final = 6


@dataclass(frozen=True)
class IndexParams:
    """Build parameters, recorded in ``meta.json``."""

    k: int = 11
    alphabet: str = "protein"
    t_base: float = 1 / 1000
    n_min: int = 8
    t_cap: float = 0.2
    oversample: float = 4.0
    mask_adapters: bool = True
    max_groups: int = 64
    fp_bits: int = 16
    t_dense: float = 0.0  # 0: no dense tier
    batch_residues: int = 20_000_000


def _mix(x: NDArray[np.uint64]) -> NDArray[np.uint64]:
    """splitmix64 finaliser, elementwise (uint64 arrays wrap)."""
    x = x ^ (x >> np.uint64(30))
    x = x * np.uint64(0xBF58476D1CE4E5B9)
    x = x ^ (x >> np.uint64(27))
    x = x * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


def _set_hashes(offsets: NDArray[np.int64], values: NDArray[np.uint64]) -> NDArray[np.uint64]:
    """Two 64-bit content hashes per (non-empty) CSR set, independent of value order."""
    lengths = np.diff(offsets).astype(np.uint64)
    out = np.empty((len(lengths), 2), dtype=np.uint64)
    for i, seed in enumerate(SET_SEEDS):
        summed = (
            np.add.reduceat(_mix(values ^ np.uint64(seed)), offsets[:-1])
            if len(lengths)
            else lengths
        )
        out[:, i] = _mix(summed + lengths * np.uint64(seed))
    return out


def _dense_ids(hashes: NDArray[np.uint64]) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Ids of distinct rows of ``hashes`` (n, 2), numbered in row order; and each id's first
    row."""
    order = np.lexsort((hashes[:, 1], hashes[:, 0]))
    ordered = hashes[order]
    new = np.ones(len(order), dtype=bool)
    new[1:] = (ordered[1:] != ordered[:-1]).any(axis=1)
    ids = np.empty(len(order), dtype=np.int64)
    ids[order] = np.cumsum(new) - 1
    return ids, order[new]


def _gather(
    offsets: NDArray[np.int64], values: NDArray[np.uint64], rows: NDArray[np.int64]
) -> tuple[NDArray[np.int64], NDArray[np.uint64]]:
    """CSR ``rows`` of (``offsets``, ``values``), as a new CSR pair."""
    lengths = offsets[rows + 1] - offsets[rows]
    out = np.concatenate([[0], np.cumsum(lengths)]).astype(np.int64)
    index = np.arange(out[-1]) - np.repeat(out[:-1] - offsets[rows], lengths)
    return out, values[index]


def packed_layout(max_hash: int, fp_bits: int, n_rows: int) -> dict[str, int]:
    """Layout of a ``PackedTable`` of ``n_rows`` (hash, value) rows, hashes <= ``max_hash``."""
    lead = 64 - max_hash.bit_length()
    fp_bits = min(fp_bits, 64 - lead)
    bucket_bits = min(max((n_rows - 1).bit_length() - KEYS_PER_BUCKET_BITS, 1), 64 - lead - fp_bits)
    return {"max_hash": max_hash, "lead_bits": lead, "bucket_bits": bucket_bits, "fp_bits": fp_bits}


def key_shift(layout: dict[str, int]) -> int:
    """Hash bits below the stored key (bucket + fingerprint bits)."""
    return 64 - layout["lead_bits"] - layout["bucket_bits"] - layout["fp_bits"]


class PackedPart(NamedTuple):
    """Keys of one hash range and their value sets, deduplicated within the range."""

    keys: NDArray[np.uint64]  # sorted, distinct
    set_ids: NDArray[np.int64]  # local set of each key
    set_offsets: NDArray[np.int64]
    set_values: NDArray[np.uint64]
    set_hashes: NDArray[np.uint64]  # (sets, 2), from _set_hashes

    @classmethod
    def build(
        cls, hashes: NDArray[np.uint64], values: NDArray[np.uint64], layout: dict[str, int]
    ) -> Self:
        """Key (hash, value) rows under ``layout``; colliding hashes pool their values."""
        pairs = (
            pl.DataFrame({"key": hashes >> np.uint64(key_shift(layout)), "value": values})
            .unique()
            .sort("key", "value")
        )
        key, value = pairs["key"].to_numpy(), pairs["value"].to_numpy()
        new = np.ones(len(key), dtype=bool)
        new[1:] = key[1:] != key[:-1]
        key_offsets = np.append(np.flatnonzero(new), len(key)).astype(np.int64)
        hashed = _set_hashes(key_offsets, value)
        set_ids, first = _dense_ids(hashed)
        set_offsets, set_values = _gather(key_offsets, value, first)
        return cls(key[new], set_ids, set_offsets, set_values, hashed[first])

    def save(self, prefix: str) -> None:
        for field in self._fields:
            np.save(f"{prefix}.{field}.npy", getattr(self, field))

    @classmethod
    def load(cls, prefix: str) -> Self:
        return cls(*(np.load(f"{prefix}.{f}.npy", mmap_mode="r") for f in cls._fields))


class SetSlice(NamedTuple):
    """Distinct value sets of one set-hash range, in hash order (``dedup_sets``), and the
    map from each part's local sets in that range to them."""

    offsets: NDArray[np.int64]
    values: NDArray[np.uint64]
    maps: NDArray[np.int64]  # every part's local sets in the range, part after part
    map_offsets: NDArray[np.int64]  # (parts + 1): each part's run in ``maps``

    def save(self, prefix: str) -> None:
        for field in self._fields:
            np.save(f"{prefix}.{field}.npy", getattr(self, field))

    @classmethod
    def load(cls, prefix: str) -> Self:
        return cls(*(np.load(f"{prefix}.{f}.npy", mmap_mode="r") for f in cls._fields))


def dedup_sets(parts: Sequence[PackedPart], lo: int, hi: int | None) -> SetSlice:
    """Distinct sets of ``parts`` whose first hash word lies in [``lo``, ``hi``) (``hi``
    None: to the end). Parts keep their sets in hash order, so each contributes one run."""
    runs = []
    for part in parts:
        first_word = part.set_hashes[:, 0]
        a = int(np.searchsorted(first_word, np.uint64(lo)))
        b = len(first_word) if hi is None else int(np.searchsorted(first_word, np.uint64(hi)))
        runs.append((part, a, b))
    hashes = np.concatenate([np.zeros((0, 2), np.uint64), *(p.set_hashes[a:b] for p, a, b in runs)])
    offsets, values, base = [], [], 0
    for part, a, b in runs:
        run = np.asarray(part.set_offsets[a : b + 1], dtype=np.int64)
        offsets.append(run[:-1] - run[0] + base)
        values.append(np.asarray(part.set_values[run[0] : run[-1]]))
        base += int(run[-1] - run[0])
    ids, first = _dense_ids(hashes)
    set_offsets, set_values = _gather(
        np.concatenate([np.zeros(0, np.int64), *offsets, [base]]).astype(np.int64),
        np.concatenate([np.zeros(0, np.uint64), *values]),
        first,
    )
    map_offsets = np.concatenate([[0], np.cumsum([b - a for _, a, b in runs])]).astype(np.int64)
    return SetSlice(set_offsets, set_values, ids, map_offsets)


@dataclass(frozen=True)
class PackedTable:
    """Hash -> value set, storing only ``fp_bits`` fingerprint bits per key.

    Hashes are at most ``max_hash``, so their top ``lead_bits`` are zero. The next
    ``bucket_bits`` address ``offsets`` directly and the next ``fp_bits`` are stored.
    Keys whose fingerprints collide share the union of their value sets. Each distinct
    value set is stored once (``set_offsets``/``set_values`` in CSR form), numbered in
    order of a 128-bit hash of its content, so a table packed in hash ranges
    (``PackedPart``) and concatenated is identical to one packed at once.
    """

    max_hash: int
    lead_bits: int
    bucket_bits: int
    fp_bits: int
    offsets: NDArray[np.unsignedinteger]
    fingerprints: NDArray[np.unsignedinteger]
    set_ids: NDArray[np.unsignedinteger]
    set_offsets: NDArray[np.unsignedinteger]
    set_values: NDArray[np.unsignedinteger]

    @property
    def shift(self) -> int:
        return 64 - self.lead_bits - self.bucket_bits - self.fp_bits

    @classmethod
    def build(
        cls, hashes: NDArray[np.uint64], values: NDArray[np.uint64], max_hash: int, fp_bits: int
    ) -> Self:
        """Pack (hash, value) pairs; every hash must be at most ``max_hash``."""
        layout = packed_layout(max_hash, fp_bits, len(hashes))
        return cls.concat(layout, [PackedPart.build(hashes, values, layout)])

    @classmethod
    def concat(
        cls,
        layout: dict[str, int],
        parts: Sequence[PackedPart],
        slices: Sequence[SetSlice] | None = None,
    ) -> Self:
        """Join parts built under ``layout`` over consecutive, disjoint key ranges.

        ``slices`` are the parts' sets deduplicated over consecutive set-hash ranges
        covering all hashes (``dedup_sets``, which can run as separate jobs); without them
        all sets are deduplicated here.
        """
        fp_bits, bucket_bits = layout["fp_bits"], layout["bucket_bits"]
        slices = slices if slices is not None else [dedup_sets(parts, 0, None)]
        set_base = np.cumsum([0] + [len(s.offsets) - 1 for s in slices])
        value_base = np.cumsum([0] + [len(s.values) for s in slices])
        set_offsets = np.concatenate(
            [
                *(
                    np.asarray(s.offsets[:-1]) + b
                    for s, b in zip(slices, value_base[:-1], strict=True)
                ),
                [value_base[-1]],
            ]
        ).astype(np.int64)
        set_values = np.concatenate([np.zeros(0, np.uint64), *(s.values for s in slices)])

        n_keys = sum(len(p.keys) for p in parts)
        fingerprints = np.empty(n_keys, dtype=_smallest(2**fp_bits - 1))
        set_ids = np.empty(n_keys, dtype=_smallest(int(set_base[-1])))
        counts = np.zeros(2**bucket_bits, dtype=np.int64)
        at, last = 0, -1
        for i, part in enumerate(parts):
            keys = np.asarray(part.keys)
            if len(keys) and int(keys[0]) <= last:
                raise ValueError("parts must cover increasing, disjoint key ranges")
            # Local sets are in hash order, so the slices' runs follow each other.
            to_global = np.concatenate(
                [
                    np.zeros(0, np.int64),
                    *(
                        b + np.asarray(s.maps[s.map_offsets[i] : s.map_offsets[i + 1]])
                        for s, b in zip(slices, set_base[:-1], strict=True)
                    ),
                ]
            )
            if len(to_global) != len(part.set_hashes):
                raise ValueError("set slices do not cover every set of every part")
            fingerprints[at : at + len(keys)] = keys & np.uint64(2**fp_bits - 1)
            set_ids[at : at + len(keys)] = to_global[np.asarray(part.set_ids)]
            buckets, n = np.unique(keys >> np.uint64(fp_bits), return_counts=True)
            counts[buckets.astype(np.int64)] += n
            at, last = at + len(keys), int(keys[-1]) if len(keys) else last
        return cls(
            **layout,
            offsets=np.concatenate([[0], np.cumsum(counts)]).astype(_smallest(n_keys)),
            fingerprints=fingerprints,
            set_ids=set_ids,
            set_offsets=set_offsets.astype(_smallest(len(set_values))),
            set_values=set_values.astype(_smallest(int(set_values.max(initial=0)))),
        )

    def lookup(self, hashes: NDArray[np.uint64]) -> NDArray[np.int64]:
        """Set id per hash, or -1 (false hits at ~keys per bucket * 2**-fp_bits per lookup)."""
        hashes = np.asarray(hashes, dtype=np.uint64)
        stored = self.stored_keys
        if not len(stored):
            return np.full(len(hashes), -1, dtype=np.int64)
        keys = hashes >> np.uint64(self.shift)
        i = np.minimum(np.searchsorted(stored, keys), len(stored) - 1)
        hit = (stored[i] == keys) & (hashes <= np.uint64(self.max_hash))
        return np.where(hit, self.set_ids[i].astype(np.int64), -1)

    def values(self, set_id: int) -> NDArray[np.uint64]:
        """Values of one set."""
        return self.set_values[self.set_offsets[set_id] : self.set_offsets[set_id + 1]]

    @cached_property
    def stored_keys(self) -> NDArray[np.uint64]:
        """Full sorted keys (bucket and fingerprint), materialised once for ``lookup``."""
        # ponytail: 8 bytes per key in RAM; the Rust lookup will walk offsets directly.
        # Shifted and filled in place: no second key-sized temporary (3.4e9 keys = 27 GB each).
        keys = np.repeat(
            np.arange(len(self.offsets) - 1, dtype=np.uint64),
            np.diff(self.offsets).astype(np.int64),
        )
        keys <<= np.uint64(self.fp_bits)
        keys |= self.fingerprints  # buffered cast, no uint64 copy of the fingerprints
        return keys

    def save(self, directory: Path, name: str) -> dict[str, int]:
        """Write one ``.npy`` per array; return the scalar layout for ``meta.json``."""
        for field in PACKED_FIELDS:
            np.save(directory / f"{name}.{field}.npy", getattr(self, field))
        return {f: getattr(self, f) for f in ("max_hash", "lead_bits", "bucket_bits", "fp_bits")}

    @classmethod
    def load(cls, directory: Path, name: str, layout: dict[str, int], *, mmap: bool = True) -> Self:
        """Memory-map a table written by ``save``, or read it into memory (``mmap=False``)."""
        return cls(
            layout["max_hash"],
            layout["lead_bits"],
            layout["bucket_bits"],
            layout["fp_bits"],
            *(
                np.load(directory / f"{name}.{f}.npy", mmap_mode="r" if mmap else None)
                for f in PACKED_FIELDS
            ),
        )

    def nbytes(self) -> int:
        return sum(getattr(self, f).nbytes for f in PACKED_FIELDS)


@dataclass(frozen=True)
class Index:
    """A built index: unit table, tiers and the metadata written by ``build_index``."""

    meta: dict[str, Any]
    units: pl.DataFrame
    tier2: PackedTable
    dense: PackedTable | None = None

    @classmethod
    def load(cls, directory: str | Path, *, mmap: bool = True) -> Self:
        """Read an index; its tiers are memory-mapped unless ``mmap`` is False."""
        directory = Path(directory)
        meta = json.loads((directory / "meta.json").read_text())
        return cls(
            meta=meta,
            units=pl.read_parquet(directory / "units.parquet"),
            tier2=PackedTable.load(directory, "tier2", meta["tier2"], mmap=mmap),
            dense=(
                PackedTable.load(directory, "dense", meta["dense"], mmap=mmap)
                if "dense" in meta
                else None
            ),
        )


def _smallest(max_value: int) -> np.dtype[np.unsignedinteger]:
    """Smallest unsigned dtype holding ``max_value`` (at least 8 bits)."""
    return np.min_scalar_type(max(max_value, 1))


def adapter_peptides() -> list[str]:
    """Stop-free ``MASK_WIDTH``-mers of the six-frame translations of ``ADAPTERS``."""
    return sorted(
        {
            frame[i : i + MASK_WIDTH]
            for adapter in ADAPTERS
            for frame in (f.decode() for f in _core.translate_frames(adapter.encode()))
            for i in range(len(frame) - MASK_WIDTH + 1)
            if "*" not in frame[i : i + MASK_WIDTH]
        }
    )


def mask_adapters(sequence: pl.Expr) -> pl.Expr:
    """Replace adapter peptides with ``X`` so no k-mer spans them."""
    return sequence.str.replace_all("|".join(adapter_peptides()), "X" * MASK_WIDTH)


def load_members(path: str | Path) -> pl.DataFrame:
    """Read a members table and number units (by ``cluster_rep``) and members from 0.

    ``counts`` marks the members ``p_in`` is computed over: full-length members, or all
    members of units that have none.
    """
    return (
        pl.read_parquet(path, columns=["protein_id", "cluster_rep", "full_length", "sequence"])
        .sort("cluster_rep", "protein_id")
        .with_columns(
            unit=(pl.col("cluster_rep").rank("dense") - 1).cast(pl.UInt32),
            member=pl.int_range(pl.len(), dtype=pl.UInt32),
            counts=pl.col("full_length") | ~pl.col("full_length").any().over("cluster_rep"),
        )
    )


def _batches(members: pl.DataFrame, batch_residues: int) -> Iterator[pl.DataFrame]:
    """Row batches of about ``batch_residues`` residues that never split a unit."""
    batch = members.with_columns(
        cum=pl.col("sequence").str.len_bytes().cast(pl.UInt64).cum_sum()
    ).select(batch=pl.col("cum").first().over("unit") // batch_residues)
    yield from members.with_columns(batch["batch"]).partition_by(
        "batch", maintain_order=True, include_key=False
    )


def _kmers(batch: pl.DataFrame, params: IndexParams, max_hash: int) -> pl.DataFrame:
    """Distinct (unit, member, hash) of a batch, with the member's ``counts`` flag."""
    out = _core.hash_proteins(
        batch["sequence"].cast(pl.Binary).to_list(),
        params.k,
        alphabet=params.alphabet,
        max_hash=max_hash,
    )
    rows = out["seq"]
    return pl.DataFrame(
        {
            "unit": batch["unit"].to_numpy()[rows],
            "member": batch["member"].to_numpy()[rows],
            "counts": batch["counts"].to_numpy()[rows],
            "hash": out["hash"],
        }
    ).unique(["member", "hash"])


def _len_cv(
    batches: list[pl.DataFrame], params: IndexParams, kept: pl.DataFrame, members: pl.DataFrame
) -> pl.DataFrame:
    """Per unit, CV over counting members of how many ``kept`` (unit, hash) each holds."""
    max_hash = int(kept["hash"].max() or 0)  # type: ignore[arg-type]
    per_member = pl.concat(
        [
            _kmers(b.filter("counts"), params, max_hash)
            .join(kept, on=["unit", "hash"], how="semi")
            .group_by("member")
            .len("kept")
            for b in batches
        ]
    )
    n = pl.col("kept").fill_null(0)
    return (
        members.filter("counts")
        .select("unit", "member")
        .join(per_member, on="member", how="left")
        .group_by("unit")
        .agg(len_cv=(n.std(ddof=0) / n.mean()).fill_nan(0.0))
    )


def _t_g(params: IndexParams) -> pl.Expr:
    """Per-unit sampling rate from ``n_members`` and ``n_kmers``: floored for non-singletons."""
    raised = pl.min_horizontal(
        pl.lit(params.t_cap), params.oversample * params.n_min / pl.col("n_kmers")
    )
    return (
        pl.when(pl.col("n_members") > 1)
        .then(pl.max_horizontal(pl.lit(params.t_base), raised))
        .otherwise(pl.lit(params.t_base))
        .clip(upper_bound=1.0)
    )


def _n_kmers(batches: list[pl.DataFrame], params: IndexParams) -> pl.DataFrame:
    """Distinct k-mers per unit (all k-mers hashed, only counts kept; units are contiguous)."""
    return pl.concat(
        [
            pl.DataFrame(
                _core.distinct_kmers(
                    b["sequence"].cast(pl.Binary).to_list(),
                    b["unit"].to_numpy(),
                    params.k,
                    alphabet=params.alphabet,
                )
            ).rename({"group": "unit"})
            for b in batches
        ]
    )


def _unit_table(members: pl.DataFrame, n_kmers: pl.DataFrame, params: IndexParams) -> pl.DataFrame:
    units = (
        members.group_by("unit")
        .agg(
            pl.col("cluster_rep").first(),
            n_members=pl.len().cast(pl.UInt32),
            n_full_length=pl.col("full_length").sum().cast(pl.UInt32),
            n_counting=pl.col("counts").sum().cast(pl.UInt32),
        )
        .join(n_kmers, on="unit", how="left")
        .with_columns(pl.col("n_kmers").fill_null(0).cast(pl.UInt32))
        .sort("unit")
        .with_columns(t_g=_t_g(params))
    )
    # Exact Rust threshold rule, evaluated once per distinct t_g.
    thresholds = units.select(pl.col("t_g").unique()).with_columns(
        max_hash_g=pl.col("t_g").map_elements(_core.max_hash, return_dtype=pl.UInt64)
    )
    return units.join(thresholds, on="t_g", maintain_order="left")


def _prepare(
    members_path: str | Path, params: IndexParams
) -> tuple[pl.DataFrame, list[pl.DataFrame], dict[str, object]]:
    """Load, mask and batch members; the returned members lose their sequences (the batches
    hold the only copy)."""
    members = load_members(members_path)
    n_residues = int(members["sequence"].str.len_bytes().cast(pl.UInt64).sum())
    n_masked = 0
    if params.mask_adapters:
        masked = members.select(mask_adapters(pl.col("sequence")))["sequence"]
        n_masked = int((masked != members["sequence"]).sum())
        members = members.with_columns(sequence=masked)
    batches = list(_batches(members, params.batch_residues))
    stats: dict[str, object] = {
        "n_proteins": members.height,
        "n_residues": n_residues,
        "n_adapter_masked_proteins": n_masked,
    }
    return members.drop("sequence"), batches, stats


def _select_postings(
    presence: pl.DataFrame, units: pl.DataFrame, params: IndexParams
) -> tuple[pl.DataFrame, pl.DataFrame, dict[str, object]]:
    """Score candidates and keep postings; return (scored, postings, stats).

    ``presence`` has ``hash``, ``unit``, ``c`` (counting members holding it) and
    ``n_groups`` (units holding it at ``t_max``) for at least every candidate (unit, hash).
    """
    scored = (
        presence.join(units.select("unit", "n_counting", "max_hash_g"), on="unit")
        .filter(pl.col("hash") <= pl.col("max_hash_g"))
        .with_columns(p_in=pl.col("c") / pl.col("n_counting"))
        .with_columns(
            # k-mers seen only in partial members get half a member's weight in the score
            score=(pl.max_horizontal("p_in", 0.5 / pl.col("n_counting")) / pl.col("n_groups")).log(
                2
            ),
            pin_q=(pl.col("p_in") * (2**PIN_BITS - 1)).round().cast(pl.UInt64),
        )
        .select("hash", "unit", "p_in", "pin_q", "n_groups", "score")
    )
    # Drop promiscuous k-mers, then keep the n_min best of each floored unit's candidates.
    floored = units.select("unit", floored=pl.col("t_g") > params.t_base)
    postings = (
        scored.filter(pl.col("n_groups") <= params.max_groups)
        .join(floored, on="unit")
        .sort(["unit", "score", "hash"], descending=[False, True, False])
        .filter(~pl.col("floored") | (pl.int_range(pl.len()).over("unit") < params.n_min))
        .drop("floored")
    )
    stats: dict[str, object] = {
        "n_clusters": units.height,
        "n_singletons": int((units["n_members"] == 1).sum()),
        "n_floored": int(floored["floored"].sum()),
        "candidates_expected": float((units["t_g"] * units["n_kmers"]).sum()),
        "candidates": scored.height,
        "promiscuous_dropped": int((scored["n_groups"] > params.max_groups).sum()),
    }
    return scored, postings, stats


def _unit_pfam(pfam_path: str | Path, members: pl.DataFrame) -> pl.DataFrame:
    """Members per (``unit``, ``pfam_accession``).

    DuckDB, as the full MGnify Pfam table exceeds Polars' 2^32-row limit on one file.
    """
    with duckdb.connect() as con:
        con.register("members", members.select("protein_id", "unit").to_arrow())
        table = con.execute(
            "SELECT unit, pfam_accession, count(DISTINCT protein_id)::UINTEGER AS n_members"
            " FROM read_parquet(?) JOIN members USING (protein_id) GROUP BY ALL",
            [str(pfam_path)],
        ).pl()
    return table.cast({"unit": members["unit"].dtype})


def build_index(
    members_path: str | Path,
    out_dir: str | Path,
    params: IndexParams | None = None,
    pfam_path: str | Path | None = None,
    postings_parquet: bool = False,
) -> dict[str, object]:
    """Build an index from a members Parquet table into ``out_dir``; return its stats."""
    params = params or IndexParams()
    members, batches, stats = _prepare(members_path, params)

    # Pass 1: distinct k-mers per unit.
    units = _unit_table(members, _n_kmers(batches, params), params)
    t_max_hash = int(units["max_hash_g"].max())  # type: ignore[arg-type]
    thresholds = units.select("unit", "max_hash_g")

    # Pass 2: candidate hashes; pass 3: presence of candidates in every unit.
    candidates = pl.concat(
        [
            _kmers(b, params, t_max_hash)
            .join(thresholds, on="unit")
            .filter(pl.col("hash") <= pl.col("max_hash_g"))
            .select("hash")
            .unique()
            for b in batches
        ]
    ).unique()
    presence = pl.concat(
        [
            _kmers(b, params, t_max_hash)
            .join(candidates, on="hash", how="semi")
            .group_by("unit", "hash")
            .agg(c=pl.col("counts").sum())
            for b in batches
        ]
    ).with_columns(n_groups=pl.len().over("hash").cast(pl.UInt32))
    _, postings, posting_stats = _select_postings(presence, units, params)
    stats |= posting_stats
    units = units.join(
        _len_cv(batches, params, postings.select("unit", "hash"), members),
        on="unit",
        how="left",
        maintain_order="left",
    )
    dense = None
    if params.t_dense > 0:
        # Pass 4: every unit's k-mers at max(t_dense, t_g), so the dense set of each unit
        # contains its tier-2 set, with p_in, minus promiscuous ones.
        units = units.with_columns(
            max_hash_dense=pl.col("max_hash_g").clip(lower_bound=_core.max_hash(params.t_dense))
        )
        dense_hash = int(units["max_hash_dense"].max())  # type: ignore[arg-type]
        dense = (
            pl.concat(
                [
                    _kmers(b, params, dense_hash)
                    .group_by("unit", "hash")
                    .agg(c=pl.col("counts").sum())
                    for b in batches
                ]
            )
            .filter(pl.len().over("hash") <= params.max_groups)
            .join(units.select("unit", "n_counting", "max_hash_dense"), on="unit")
            .filter(pl.col("hash") <= pl.col("max_hash_dense"))
            .with_columns(p_in=pl.col("c") / pl.col("n_counting"))
            .select(
                "hash",
                "unit",
                "p_in",
                pin_q=(pl.col("p_in") * (2**PIN_BITS - 1)).round().cast(pl.UInt64),
            )
        )
        units = units.join(
            _len_cv(batches, params, dense.select("unit", "hash"), members).rename(
                {"len_cv": "len_cv_dense"}
            ),
            on="unit",
            how="left",
            maintain_order="left",
        )
    return finish_index(
        out_dir,
        params,
        units,
        postings,
        stats,
        unit_pfam=None if pfam_path is None else _unit_pfam(pfam_path, members),
        dense=dense,
        postings_parquet=postings_parquet,
    )


def finish_index(
    out_dir: str | Path,
    params: IndexParams,
    units: pl.DataFrame,
    postings: pl.DataFrame,
    stats: dict[str, object],
    unit_pfam: pl.DataFrame | None = None,
    dense: pl.DataFrame | None = None,
    postings_parquet: bool = False,
) -> dict[str, object]:
    """Drop units without a posting (they can never be hit), renumber the rest in ``unit``
    order and write the index."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ids = (
        postings.select("unit")
        .unique()
        .sort("unit")
        .with_columns(new=pl.int_range(pl.len(), dtype=pl.UInt32))
    )

    def renumber(table: pl.DataFrame) -> pl.DataFrame:
        return (
            table.join(ids, on="unit", maintain_order="left")
            .with_columns(unit=pl.col("new"))
            .drop("new")
        )

    if unit_pfam is not None:
        renumber(unit_pfam).sort("unit", "pfam_accession").write_parquet(out / "unit_pfam.parquet")
    units = renumber(
        units.select(
            "unit",
            "cluster_rep",
            "n_members",
            "n_kmers",
            "t_g",
            "max_hash_g",
            "^len_cv.*$",
            *(("max_hash_dense",) if dense is not None else ()),
        )
    )
    return write_index(
        out,
        params,
        units,
        renumber(postings),
        stats,
        dense=None if dense is None else renumber(dense),
        postings_parquet=postings_parquet,
    )


def _pin_hist(postings: pl.DataFrame, n_units: int) -> pl.Series:
    """Per unit, how many kept k-mers sit at each quantised ``p_in`` level."""
    units = postings["unit"].to_numpy()
    hist = np.zeros((n_units, 2**PIN_BITS), dtype=_smallest(int(np.bincount(units).max(initial=0))))
    np.add.at(hist, (units, postings["pin_q"].to_numpy()), 1)
    return pl.Series(hist)


def _pin_sum(postings: pl.DataFrame, n_units: int) -> pl.Series:
    """Per unit, the sum of ``p_in`` over kept k-mers: the kept k-mers of an average member."""
    return pl.Series(
        np.bincount(postings["unit"].to_numpy(), postings["p_in"].to_numpy(), minlength=n_units)
    ).cast(pl.Float32)


def unit_columns(units: pl.DataFrame, postings: pl.DataFrame) -> pl.DataFrame:
    """Add ``pin_hist``, ``pin_sum``, ``m_g`` and a filled f32 ``len_cv`` to ``units``
    (``unit`` 0..n-1) from their ``postings``."""
    units = units.with_columns(
        pin_hist=_pin_hist(postings, units.height), pin_sum=_pin_sum(postings, units.height)
    )
    # Imported sketches have no members, so no length spread; units without members get 0.
    units = units.with_columns(
        (pl.col("len_cv") if "len_cv" in units.columns else pl.lit(0.0))
        .fill_null(0.0)
        .cast(pl.Float32)
        .alias("len_cv")
    )
    return units.join(
        postings.group_by("unit").agg(m_g=pl.len().cast(pl.UInt32)),
        on="unit",
        how="left",
        maintain_order="left",
    ).with_columns(pl.col("m_g").fill_null(0))


def write_meta(
    out: Path,
    params: IndexParams,
    tables: dict[str, PackedTable],
    stats: dict[str, object],
    hash_scheme: str = "kfp",
) -> dict[str, object]:
    """Save the tables and ``meta.json``; ``stats`` gains the tier-2 sizes. Returns it."""
    tier2 = tables["tier2"]
    stats = {
        **stats,
        "tier2_keys": len(tier2.fingerprints),
        "tier2_sets": len(tier2.set_offsets) - 1,
        "tier2_bytes": tier2.nbytes(),
        "tier2_bytes_per_hash": tier2.nbytes() / max(stats["distinct_hashes"], 1),  # type: ignore[call-overload]
    }
    meta = {
        "format": 2,
        "hash": hash_scheme,
        "params": asdict(params),
        **{name: table.save(out, name) for name, table in tables.items()},
        "stats": stats,
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    return stats


def write_index(
    out: Path,
    params: IndexParams,
    units: pl.DataFrame,
    postings: pl.DataFrame,
    stats: dict[str, object],
    hash_scheme: str = "kfp",
    dense: pl.DataFrame | None = None,
    postings_parquet: bool = False,
) -> dict[str, object]:
    """Pack the tiers and write the index files.

    ``units`` needs ``unit`` (0..n-1), ``t_g`` and ``max_hash_g``; ``postings`` holds one row
    per kept (hash, unit), each with ``hash <= max_hash_g``: ``p_in``, ``pin_q``,
    ``n_groups`` and ``score``. ``hash_scheme`` tells the query how to hash reads. ``dense``
    (``hash``, ``unit``, ``p_in``, ``pin_q``) becomes the dense table. Returns ``stats``
    extended with sizes. ``postings_parquet`` also writes ``postings`` for inspection.
    """
    t_max_hash = int(units["max_hash_g"].max())  # type: ignore[arg-type]
    postings = postings.sort("hash", "unit")
    n_units = units.height
    units = unit_columns(units, postings)

    tier2 = PackedTable.build(
        postings["hash"].to_numpy(),
        (postings["unit"].cast(pl.UInt64).to_numpy() << np.uint64(PIN_BITS))
        | postings["pin_q"].to_numpy(),
        t_max_hash,
        params.fp_bits,
    )
    tables = {"tier2": tier2}
    if dense is not None:
        tables["dense"] = PackedTable.build(
            dense["hash"].to_numpy(),
            (dense["unit"].cast(pl.UInt64).to_numpy() << np.uint64(PIN_BITS))
            | dense["pin_q"].to_numpy(),
            int(units["max_hash_dense"].max()),  # type: ignore[arg-type]
            params.fp_bits,
        )
        units = units.join(
            dense.group_by("unit").agg(m_dense=pl.len().cast(pl.UInt32)),
            on="unit",
            how="left",
            maintain_order="left",
        ).with_columns(
            pl.col("m_dense").fill_null(0),
            pl.col("len_cv_dense").fill_null(0.0).cast(pl.Float32),
            pin_hist_dense=_pin_hist(dense, units.height),
            pin_sum_dense=_pin_sum(dense, units.height),
        )
        stats["dense_postings"] = dense.height
        stats["dense_bytes"] = tables["dense"].nbytes()

    units.write_parquet(out / "units.parquet")
    if postings_parquet:
        postings.write_parquet(out / "postings.parquet")
    stats = {
        **stats,
        "n_units": n_units,
        "t_max": float(units["t_g"].max()),  # type: ignore[arg-type]
        "postings": postings.height,
        "distinct_hashes": postings["hash"].n_unique(),
    }
    return write_meta(out, params, tables, stats, hash_scheme)
