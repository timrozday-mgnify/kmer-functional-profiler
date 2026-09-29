from collections.abc import Sequence
from os import PathLike
from typing import Self

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
