from collections.abc import Sequence
from os import PathLike
from typing import Protocol, Self

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
    seqs: Sequence[bytes], k: int, *, alphabet: str = "protein", max_hash: int = ...
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
def gather(
    units: NDArray[np.uint32], hashes: NDArray[np.uint64], t_g: NDArray[np.float64]
) -> dict[str, NDArray[np.uint32]]: ...

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
    ) -> None: ...
    def __iter__(self) -> Self: ...
    def __next__(self) -> dict[str, NDArray[np.uint64] | NDArray[np.uint8]]: ...
    @property
    def n_reads(self) -> int: ...
