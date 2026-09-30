"""FASTQ streaming on the fixtures in tests/data (see scripts/make_fixtures.py)."""

import csv
import gzip
from pathlib import Path

import numpy as np
import pytest

from kmer_functional_profiler import _core, reference

DATA = Path(__file__).resolve().parent.parent / "data"
R1, R2 = DATA / "reads_1.fastq.gz", DATA / "reads_2.fastq.gz"


def read_fastq(path: Path) -> list[bytes]:
    with gzip.open(path, "rb") as f:
        return f.read().splitlines()[1::4]


def read_proteins() -> dict[str, str]:
    lines = (DATA / "proteins.faa").read_text().splitlines()
    return {name[1:]: seq for name, seq in zip(lines[::2], lines[1::2], strict=True)}


def concat(batches: list[dict[str, np.ndarray]]) -> dict[str, list[int]]:
    return {name: np.concatenate([b[name] for b in batches]).tolist() for name in batches[0]}


@pytest.mark.parametrize("batch_reads", [1, 7, 100_000])
def test_paired_hits_match_reference(batch_reads: int) -> None:
    mate1, mate2 = read_fastq(R1), read_fastq(R2)
    hits = _core.FastxHits(R1, R2, k=8, frames="all", batch_reads=batch_reads)
    got = concat(list(hits))
    seqs = [s for pair in zip(mate1, mate2, strict=True) for s in pair]
    reads = [i for i in range(len(mate1)) for _ in range(2)]
    expected = reference.hash_dna(seqs, 8, frames="all", reads=reads, mates=[0, 1] * len(mate1))
    assert got == expected
    assert hits.n_reads == len(mate1)


@pytest.mark.parametrize("batch_reads", [5000, 100_000])
def test_threaded_batches_keep_input_order(tmp_path: Path, batch_reads: int) -> None:
    # Enough pairs that a batch is split over threads (4096 reads per thread at least).
    rng = np.random.default_rng(0)
    n = 3 * 4096 + 5
    seqs = [
        bytes(s) for s in np.frombuffer(b"ACGT", dtype=np.uint8)[rng.integers(0, 4, (2 * n, 140))]
    ]
    for mate in (0, 1):
        with gzip.open(tmp_path / f"r{mate}.fastq.gz", "wb") as f:
            f.writelines(b"@r\n%s\n+\n%s\n" % (s, b"I" * 140) for s in seqs[mate::2])
    hits = _core.FastxHits(
        tmp_path / "r0.fastq.gz", tmp_path / "r1.fastq.gz", k=8, batch_reads=batch_reads
    )
    got = concat(list(hits))
    expected = {name: col.tolist() for name, col in _core.hash_dna(seqs, 8).items()}
    expected["mate"] = [r % 2 for r in expected["read"]]
    expected["read"] = [r // 2 for r in expected["read"]]
    assert got == expected
    assert hits.n_reads == n


def test_true_frames_translate_to_the_source_protein() -> None:
    reads = {0: read_fastq(R1), 1: read_fastq(R2)}
    proteins = read_proteins()
    with (DATA / "truth.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert {int(r["frame"]) for r in rows} == set(range(6))
    checked = 0
    for row in rows:
        read = reads[int(row["mate"])][int(row["read"])]
        translation = _core.translate_frames(read)[int(row["frame"])].decode()
        start, end = int(row["aa_start"]), int(row["aa_end"])
        source = proteins[row["protein"]]
        if row["case"] in ("clean", "synonymous"):
            assert translation == source[start:end]
            checked += 1
        elif row["case"] == "stop":
            assert translation.startswith(source[start:] + "*")
            stopfree = _core.hash_dna([read], 8)["frame"]
            assert int(row["frame"]) not in stopfree.tolist()
            checked += 1
    assert checked > 40


def test_mate_count_mismatch_is_an_error(tmp_path: Path) -> None:
    short = tmp_path / "short.fastq"
    short.write_text("@r\nACGTACGTACGT\n+\nIIIIIIIIIIII\n")
    with pytest.raises(ValueError, match="mate files"):
        list(_core.FastxHits(R1, short, k=5))


def test_missing_file_is_an_os_error(tmp_path: Path) -> None:
    with pytest.raises(OSError):
        _core.FastxHits(tmp_path / "missing.fastq", k=5)


def test_invalid_parameters_are_value_errors() -> None:
    with pytest.raises(ValueError, match="k = 13"):
        _core.hash_dna([b"ACGT"], 13)
    with pytest.raises(ValueError, match="alphabet"):
        _core.hash_dna([b"ACGT"], 5, alphabet="nope")
    with pytest.raises(ValueError, match="genetic code"):
        _core.hash_dna([b"ACGT"], 5, genetic_code=2)
