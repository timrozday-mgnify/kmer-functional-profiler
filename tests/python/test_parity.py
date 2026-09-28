"""Rust kernels agree with the pure-Python reference twins on random inputs."""

import math
import random

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from kmer_functional_profiler import _core, reference

ALPHABETS = list(reference.BITS)
dna = st.text(alphabet="ACGTNacgtu", max_size=240).map(str.encode)
protein = st.text(alphabet="ACDEFGHIKLMNPQRSTVWYXBZ*acdefghik", max_size=120).map(str.encode)
thresholds = st.one_of(st.just(reference.MAX_HASH), st.integers(0, reference.MAX_HASH))


@st.composite
def alphabet_and_k(draw: st.DrawFn) -> tuple[str, int]:
    alphabet = draw(st.sampled_from(ALPHABETS))
    return alphabet, draw(st.integers(1, 64 // reference.BITS[alphabet]))


def as_lists(columns: dict[str, np.ndarray]) -> dict[str, list[int]]:
    return {name: col.tolist() for name, col in columns.items()}


@given(dna, st.sampled_from([11, 4]))
def test_translate_frames(seq: bytes, code: int) -> None:
    assert _core.translate_frames(seq, code) == reference.six_frames(seq, code)


@settings(max_examples=200)
@given(
    st.lists(dna, max_size=4),
    alphabet_and_k(),
    st.sampled_from([11, 4]),
    st.sampled_from(["stopfree", "all"]),
    thresholds,
)
def test_hash_dna(
    seqs: list[bytes], ak: tuple[str, int], code: int, frames: str, threshold: int
) -> None:
    alphabet, k = ak
    kwargs = {
        "alphabet": alphabet,
        "genetic_code": code,
        "frames": frames,
        "max_hash": threshold,
    }
    got = as_lists(_core.hash_dna(seqs, k, **kwargs))
    assert got == reference.hash_dna(seqs, k, **kwargs)


@given(st.lists(protein, max_size=4), alphabet_and_k(), thresholds)
def test_hash_proteins(seqs: list[bytes], ak: tuple[str, int], threshold: int) -> None:
    alphabet, k = ak
    got = as_lists(_core.hash_proteins(seqs, k, alphabet=alphabet, max_hash=threshold))
    assert got == reference.hash_proteins(seqs, k, alphabet=alphabet, max_hash=threshold)


@given(st.floats(allow_nan=True, allow_infinity=True))
def test_max_hash(fraction: float) -> None:
    assert _core.max_hash(fraction) == reference.max_hash(fraction)


def test_max_hash_keeps_the_requested_fraction() -> None:
    threshold = _core.max_hash(0.01)
    seq = "".join(random.Random(0).choices("ACDEFGHIKLMNPQRSTVWY", k=100_000)).encode()
    hashes = _core.hash_proteins([seq], 8)["hash"]
    kept = np.mean(hashes <= np.uint64(threshold))
    assert math.isclose(kept, 0.01, abs_tol=0.002)
