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
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from functools import cached_property
from pathlib import Path
from typing import Any, Final, Self

import numpy as np
import polars as pl
from numpy.typing import NDArray

from kmer_functional_profiler import _core

U64_MAX: Final = 2**64 - 1
PIN_BITS: Final = 4  # tier-2 values are unit << PIN_BITS | quantised p_in
PACKED_FIELDS: Final = ("offsets", "fingerprints", "set_ids", "set_offsets", "set_values")
KEYS_PER_BUCKET_BITS: Final = 2  # ~4 keys per offsets bucket
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


@dataclass(frozen=True)
class PackedTable:
    """Hash -> value set, storing only ``fp_bits`` fingerprint bits per key.

    Hashes are at most ``max_hash``, so their top ``lead_bits`` are zero. The next
    ``bucket_bits`` address ``offsets`` directly and the next ``fp_bits`` are stored.
    Keys whose fingerprints collide share the union of their value sets. Each distinct
    value set is stored once (``set_offsets``/``set_values`` in CSR form).
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
        lead = 64 - max_hash.bit_length()
        fp_bits = min(fp_bits, 64 - lead)
        n_keys = len(np.unique(hashes))
        bucket_bits = min(
            max((n_keys - 1).bit_length() - KEYS_PER_BUCKET_BITS, 1), 64 - lead - fp_bits
        )
        shift = 64 - lead - bucket_bits - fp_bits
        keyed = (
            pl.DataFrame({"key": hashes >> np.uint64(shift), "value": values})
            .unique()
            .group_by("key")
            .agg(pl.col("value").sort())
            .sort("key")
            # ponytail: string-keyed set dedup; hash the lists if this is slow at scale
            .with_columns(
                name=pl.col("value").list.eval(pl.element().cast(pl.String)).list.join(",")
            )
            .with_columns(set_id=(pl.col("name").rank("dense") - 1).cast(pl.UInt32))
        )
        sets = keyed.unique("set_id").sort("set_id")
        keys = keyed["key"].to_numpy()
        lengths = sets["value"].list.len().to_numpy()
        set_values = sets["value"].explode(empty_as_null=False).to_numpy()
        return cls(
            max_hash=max_hash,
            lead_bits=lead,
            bucket_bits=bucket_bits,
            fp_bits=fp_bits,
            offsets=np.searchsorted(
                keys >> np.uint64(fp_bits), np.arange(2**bucket_bits + 1, dtype=np.uint64)
            ).astype(_smallest(len(keys))),
            fingerprints=(keys & np.uint64(2**fp_bits - 1)).astype(_smallest(2**fp_bits - 1)),
            set_ids=keyed["set_id"].to_numpy().astype(_smallest(sets.height)),
            set_offsets=np.concatenate([[0], np.cumsum(lengths)]).astype(
                _smallest(len(set_values))
            ),
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
        # ponytail: 8 bytes per key in RAM; the Rust lookup will walk offsets directly
        buckets = np.repeat(
            np.arange(len(self.offsets) - 1, dtype=np.uint64),
            np.diff(self.offsets).astype(np.int64),
        )
        return (buckets << np.uint64(self.fp_bits)) | self.fingerprints.astype(np.uint64)

    def save(self, directory: Path, name: str) -> dict[str, int]:
        """Write one ``.npy`` per array; return the scalar layout for ``meta.json``."""
        for field in PACKED_FIELDS:
            np.save(directory / f"{name}.{field}.npy", getattr(self, field))
        return {f: getattr(self, f) for f in ("max_hash", "lead_bits", "bucket_bits", "fp_bits")}

    @classmethod
    def load(cls, directory: Path, name: str, layout: dict[str, int]) -> Self:
        """Memory-map a table written by ``save``."""
        return cls(
            layout["max_hash"],
            layout["lead_bits"],
            layout["bucket_bits"],
            layout["fp_bits"],
            *(np.load(directory / f"{name}.{f}.npy", mmap_mode="r") for f in PACKED_FIELDS),
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
    def load(cls, directory: str | Path) -> Self:
        directory = Path(directory)
        meta = json.loads((directory / "meta.json").read_text())
        return cls(
            meta=meta,
            units=pl.read_parquet(directory / "units.parquet"),
            tier2=PackedTable.load(directory, "tier2", meta["tier2"]),
            dense=PackedTable.load(directory, "dense", meta["dense"]) if "dense" in meta else None,
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


def build_index(
    members_path: str | Path,
    out_dir: str | Path,
    params: IndexParams | None = None,
    pfam_path: str | Path | None = None,
    postings_parquet: bool = False,
) -> dict[str, object]:
    """Build an index from a members Parquet table into ``out_dir``; return its stats."""
    params = params or IndexParams()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    members = load_members(members_path)
    n_residues = int(members["sequence"].str.len_bytes().cast(pl.UInt64).sum())
    n_masked = 0
    if params.mask_adapters:
        masked = members.select(mask_adapters(pl.col("sequence")))["sequence"]
        n_masked = int((masked != members["sequence"]).sum())
        members = members.with_columns(sequence=masked)
    batches = list(_batches(members, params.batch_residues))
    members = members.drop("sequence")  # the batches hold the only copy

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
    presence = (
        pl.concat(
            [
                _kmers(b, params, t_max_hash)
                .join(candidates, on="hash", how="semi")
                .group_by("unit", "hash")
                .agg(c=pl.col("counts").sum())
                for b in batches
            ]
        )
        .with_columns(n_groups=pl.len().over("hash").cast(pl.UInt32))
        .join(units.select("unit", "n_counting", "max_hash_g"), on="unit")
    )
    scored = (
        presence.filter(pl.col("hash") <= pl.col("max_hash_g"))
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
    units = units.join(
        _len_cv(batches, params, postings.select("unit", "hash"), members),
        on="unit",
        how="left",
        maintain_order="left",
    )
    stats: dict[str, object] = {
        "n_proteins": members.height,
        "n_clusters": units.height,
        "n_residues": n_residues,
        "n_adapter_masked_proteins": n_masked,
        "n_singletons": int((units["n_members"] == 1).sum()),
        "n_floored": floored["floored"].sum(),
        "candidates_expected": float((units["t_g"] * units["n_kmers"]).sum()),
        "candidates": scored.height,
        "promiscuous_dropped": int((scored["n_groups"] > params.max_groups).sum()),
    }
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
    # Units without a posting can never be hit: drop them and renumber the rest.
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

    if pfam_path is not None:
        renumber(
            pl.read_parquet(pfam_path, columns=["protein_id", "pfam_accession"])
            .join(members.select("protein_id", "unit"), on="protein_id")
            .group_by("unit", "pfam_accession")
            .agg(n_members=pl.col("protein_id").n_unique().cast(pl.UInt32))
        ).sort("unit", "pfam_accession").write_parquet(out / "unit_pfam.parquet")
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
    units = units.with_columns(
        pin_hist=_pin_hist(postings, units.height), pin_sum=_pin_sum(postings, units.height)
    )
    # Imported sketches have no members, so no length spread; units without members get 0.
    units = units.with_columns(
        (pl.col(c) if c in units.columns else pl.lit(0.0)).fill_null(0.0).cast(pl.Float32).alias(c)
        for c in ("len_cv", *(("len_cv_dense",) if dense is not None else ()))
    )
    units = units.join(
        postings.group_by("unit").agg(m_g=pl.len().cast(pl.UInt32)),
        on="unit",
        how="left",
        maintain_order="left",
    ).with_columns(pl.col("m_g").fill_null(0))

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
            pin_hist_dense=_pin_hist(dense, units.height),
            pin_sum_dense=_pin_sum(dense, units.height),
        )
        stats["dense_postings"] = dense.height
        stats["dense_bytes"] = tables["dense"].nbytes()

    units.write_parquet(out / "units.parquet")
    if postings_parquet:
        postings.write_parquet(out / "postings.parquet")
    distinct_hashes = postings["hash"].n_unique()
    stats = {
        **stats,
        "n_units": n_units,
        "t_max": float(units["t_g"].max()),  # type: ignore[arg-type]
        "postings": postings.height,
        "distinct_hashes": distinct_hashes,
        "tier2_keys": len(tier2.fingerprints),
        "tier2_sets": len(tier2.set_offsets) - 1,
        "tier2_bytes": tier2.nbytes(),
        "tier2_bytes_per_hash": tier2.nbytes() / max(distinct_hashes, 1),
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
