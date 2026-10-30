from collections.abc import Sequence
from os import PathLike
from typing import Any, Protocol, Self

import numpy as np
from numpy.typing import NDArray

__version__: str

def translate_frames(seq: bytes, genetic_code: int = 11) -> list[bytes]: ...
def max_hash(fraction: float) -> int: ...
def hash_dna(
    seqs: Sequence[bytes],
    k: int,
    *,
    alphabet: str = "protein",
    genetic_code: int = 11,
    frames: str = "stopfree",
    max_hash: int = ...,
) -> dict[str, NDArray[np.uint64] | NDArray[np.uint8]]: ...
def hash_proteins(
    seqs: Sequence[bytes],
    k: int,
    *,
    alphabet: str = "protein",
    max_hash: int = ...,
    min_hash: int = 0,
) -> dict[str, NDArray[np.uint64]]: ...
def distinct_kmers(
    seqs: Sequence[bytes],
    groups: NDArray[np.uint32],
    k: int,
    *,
    alphabet: str = "protein",
    max_hash: int = ...,
) -> dict[str, NDArray[np.uint32]]: ...

BLOOM_BLOCK_BYTES: int

def bloom_insert(bits: NDArray[np.uint8], hashes: NDArray[np.uint64]) -> None: ...
def bloom_contains(bits: NDArray[np.uint8], hashes: NDArray[np.uint64]) -> NDArray[np.bool_]: ...

class _Packed(Protocol):  # index.PackedTable
    @property
    def max_hash(self) -> int: ...
    @property
    def shift(self) -> int: ...
    @property
    def fp_bits(self) -> int: ...
    @property
    def offsets(self) -> NDArray[np.unsignedinteger]: ...
    @property
    def fingerprints(self) -> NDArray[np.unsignedinteger]: ...
    @property
    def set_ids(self) -> NDArray[np.unsignedinteger]: ...
    @property
    def set_offsets(self) -> NDArray[np.unsignedinteger]: ...
    @property
    def set_values(self) -> NDArray[np.unsignedinteger]: ...

def packed_lookup(table: _Packed, hashes: NDArray[np.uint64]) -> NDArray[np.int64]: ...
def unit_hits(
    table: _Packed,
    max_hash_g: NDArray[np.uint64],
    hashes: NDArray[np.uint64],
    reads: NDArray[np.uint64],
) -> dict[str, NDArray[np.uint64] | NDArray[np.uint32] | NDArray[np.uint8]]: ...
def route_hits(
    table: _Packed, reads: NDArray[np.uint64], hashes: NDArray[np.uint64]
) -> dict[str, NDArray[np.uint64] | NDArray[np.uint32]]: ...
def gather(
    units: NDArray[np.uint32], hashes: NDArray[np.uint64], t_g: NDArray[np.float64]
) -> dict[str, NDArray[np.uint32] | NDArray[np.bool_]]: ...
def components(
    units: NDArray[np.uint32], hashes: NDArray[np.uint64], n: int
) -> NDArray[np.uint32]: ...
def overlap_pairs(
    offsets: NDArray[np.uint64],
    hashes: NDArray[np.uint64],
    max_hash: NDArray[np.uint64],
    *,
    min_shared: int = 1,
    storage: bool = True,
    max_holders: int = ...,
) -> dict[str, Any]: ...
def pair_weights(
    offsets: NDArray[np.uint64],
    hashes: NDArray[np.uint64],
    max_hash: NDArray[np.uint64],
    a: NDArray[np.uint32],
    b: NDArray[np.uint32],
) -> dict[str, NDArray[np.uint32]]: ...
def kruskal(
    n: int, u: NDArray[np.uint32], v: NDArray[np.uint32], w: NDArray[np.uint32]
) -> dict[str, NDArray[np.uint32] | NDArray[np.bool_]]: ...
def tree_layout(
    n: int, left: NDArray[np.uint32], right: NDArray[np.uint32], tops: NDArray[np.uint32]
) -> dict[str, NDArray[np.uint32]]: ...
def tree_schemes(
    n: int,
    left: NDArray[np.uint32],
    right: NDArray[np.uint32],
    tops: NDArray[np.uint32],
    offsets: NDArray[np.uint64],
    hashes: NDArray[np.uint64],
    leaf_level: NDArray[np.uint8],
    level_hash: NDArray[np.uint64],
    *,
    emit: bool = False,
) -> dict[str, NDArray[np.uint64] | NDArray[np.uint32] | NDArray[np.bool_]]: ...
def overlap_components(
    offsets: NDArray[np.uint64], hashes: NDArray[np.uint64], max_hash: NDArray[np.uint64]
) -> NDArray[np.uint32]: ...
def contained_sums(
    lo: NDArray[np.uint32],
    hi: NDArray[np.uint32],
    w: NDArray[np.int64],
    q_lo: NDArray[np.uint32],
    q_hi: NDArray[np.uint32],
) -> NDArray[np.int64]: ...
def align_identity(
    a: Sequence[bytes], b: Sequence[bytes]
) -> dict[str, NDArray[np.float32] | NDArray[np.uint32]]: ...

class DistinctHashes:
    def __init__(self, max_hash: int) -> None: ...
    def add(self, hashes: NDArray[np.uint64]) -> None: ...
    def __len__(self) -> int: ...
    def take(self, i: int) -> NDArray[np.uint64]: ...

class FastxHits:
    def __init__(
        self,
        r1: str | PathLike[str],
        r2: str | PathLike[str] | None = None,
        *,
        k: int,
        alphabet: str = "protein",
        genetic_code: int = 11,
        frames: str = "stopfree",
        max_hash: int = ...,
        batch_reads: int = 100_000,
        min_qual: int = 0,
    ) -> None: ...
    def __iter__(self) -> Self: ...
    def __next__(self) -> dict[str, NDArray[np.uint64] | NDArray[np.uint8]]: ...
    @property
    def n_reads(self) -> int: ...
    @property
    def base_stats(self) -> dict[str, float]: ...

class DiamondReader:
    def __init__(
        self,
        names: Sequence[bytes],
        seqs: Sequence[bytes],
        units: NDArray[np.uint32],
        *,
        top: float = ...,
        max_bits_below: float = ...,
        unit_bits_below: float = ...,
        min_identity: float = 0.0,
    ) -> None: ...
    def parse(self, buf: bytes, batch_reads: int, eof: bool) -> tuple[dict[str, Any], int]: ...
