"""Slow pure-Python twins of the Rust kernels, used only as test oracles.

They follow the same conventions as ``crates/core`` but are written independently
(no rolling hash, no shared lookup tables) so a bug in one is unlikely to hide in both.
"""

from collections.abc import Iterable, Sequence
from typing import Final

CODE_11: Final = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE_4: Final = CODE_11[:14] + "W" + CODE_11[15:]
CODES: Final = {11: CODE_11, 4: CODE_4}
BASES: Final = "TCAG"
ALPHABETS: Final = {
    "protein": list("ACDEFGHIKLMNPQRSTVWY"),
    "murphy10": ["LVIM", "C", "A", "G", "ST", "P", "FYW", "EDNQ", "KR", "H"],
    "dayhoff": ["AGPST", "C", "DENQ", "HKR", "ILMV", "FWY"],
}
BITS: Final = {"protein": 5, "murphy10": 4, "dayhoff": 3}
MAX_HASH: Final = 2**64 - 1

_COMPLEMENT: Final = dict(zip("ACGTU", "TGCAA", strict=True))


def translate(dna: bytes, genetic_code: int = 11) -> bytes:
    """Translate every full codon; stops are ``*``, codons with non-ACGTU bases ``X``."""
    table = CODES[genetic_code]
    seq = dna.decode().upper().replace("U", "T")
    out = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if all(b in BASES for b in codon):
            a, b, c = (BASES.index(x) for x in codon)
            out.append(table[16 * a + 4 * b + c])
        else:
            out.append("X")
    return "".join(out).encode()


def reverse_complement(dna: bytes) -> bytes:
    """Reverse complement, uppercased; non-ACGT bases become ``N``."""
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed(dna.decode().upper())).encode()


def six_frames(dna: bytes, genetic_code: int = 11) -> list[bytes]:
    """Translate frames 0-2 (forward offsets) and 3-5 (reverse-complement offsets)."""
    rc = reverse_complement(dna)
    return [translate(strand[offset:], genetic_code) for strand in (dna, rc) for offset in range(3)]


def hash_kmer(kmer: int) -> int:
    """Splitmix64 finalizer."""
    z = (kmer + 0x9E3779B97F4A7C15) & MAX_HASH
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MAX_HASH
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MAX_HASH
    return z ^ (z >> 31)


def max_hash(fraction: float) -> int:
    """Threshold keeping ``fraction`` of hashes (keep iff ``hash <= max_hash``)."""
    if fraction >= 1:
        return MAX_HASH
    if not fraction > 0:  # also catches NaN
        return 0
    return max(int(fraction * 2.0**64) - 1, 0)


def _class_of(alphabet: str) -> dict[str, int]:
    return {aa: i for i, group in enumerate(ALPHABETS[alphabet]) for aa in group + group.lower()}


def protein_kmers(
    protein: bytes, k: int, alphabet: str = "protein", threshold: int = MAX_HASH
) -> list[int]:
    """Hashes of all k-mers within runs of alphabet residues, in order, that pass the threshold."""
    classes = _class_of(alphabet)
    bits = BITS[alphabet]
    hashes = []
    seq = protein.decode()
    for start in range(len(seq) - k + 1):
        window = seq[start : start + k]
        if all(aa in classes for aa in window):
            packed = 0
            for aa in window:
                packed = (packed << bits) | classes[aa]
            h = hash_kmer(packed)
            if h <= threshold:
                hashes.append(h)
    return hashes


def _segments(aa: bytes, frames: str) -> list[bytes]:
    """The parts of one translated frame that a frame mode hashes."""
    if frames == "all" or b"*" not in aa:
        return [aa]
    if frames == "stopfree":
        return []
    min_len = int(frames.partition(":")[2] or 30)  # "edges" or "edges:M"
    parts = aa.split(b"*")
    return [seg for seg in (parts[0], parts[-1]) if len(seg) >= min_len]


def mask_quality(seq: bytes, qual: bytes, min_qual: int) -> bytes:
    """Bases with Phred quality below ``min_qual`` become ``N`` (Phred+33)."""
    return bytes(b if q - 33 >= min_qual else ord("N") for b, q in zip(seq, qual, strict=True))


def hash_dna(
    seqs: Iterable[bytes],
    k: int,
    *,
    alphabet: str = "protein",
    genetic_code: int = 11,
    frames: str = "stopfree",
    max_hash: int = MAX_HASH,
    mates: Sequence[int] | None = None,
    reads: Sequence[int] | None = None,
) -> dict[str, list[int]]:
    """Twin of ``_core.hash_dna``; ``reads``/``mates`` override the per-sequence labels."""
    out: dict[str, list[int]] = {"read": [], "mate": [], "frame": [], "hash": []}
    for i, seq in enumerate(seqs):
        for frame, aa in enumerate(six_frames(seq, genetic_code)):
            for seg in _segments(aa, frames):
                for h in protein_kmers(seg, k, alphabet, max_hash):
                    out["read"].append(reads[i] if reads else i)
                    out["mate"].append(mates[i] if mates else 0)
                    out["frame"].append(frame)
                    out["hash"].append(h)
    return out


def hash_proteins(
    seqs: Iterable[bytes],
    k: int,
    *,
    alphabet: str = "protein",
    max_hash: int = MAX_HASH,
) -> dict[str, list[int]]:
    """Twin of ``_core.hash_proteins``."""
    out: dict[str, list[int]] = {"seq": [], "hash": []}
    for i, seq in enumerate(seqs):
        hashes = protein_kmers(seq, k, alphabet, max_hash)
        out["seq"] += [i] * len(hashes)
        out["hash"] += hashes
    return out
