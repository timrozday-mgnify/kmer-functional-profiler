"""Nucleotide indexes (plan, phase 14, steps 2-3): the query follows the index's alphabet."""

import json
import random
from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal
from typer.testing import CliRunner

from kmer_functional_profiler import reference
from kmer_functional_profiler.cli import app
from kmer_functional_profiler.index import DNA_K, Index, IndexParams, build_index
from kmer_functional_profiler.query import profile

RNG = random.Random(3)
GENES = ["".join(RNG.choices("ACGT", k=900)) for _ in range(3)]


def write_members(path: Path, seqs: list[str]) -> Path:
    rows = [(i, i, True, s) for i, s in enumerate(seqs)]
    pl.DataFrame(
        rows, schema=["protein_id", "cluster_rep", "full_length", "sequence"], orient="row"
    ).write_parquet(path)
    return path


def write_reads(path: Path, gene: str, n: int, strand: str = "both", length: int = 150) -> Path:
    """Reads from random starts on ``gene`` (the same starts for every ``strand``): forward,
    reverse complemented, or every other one reversed (``both``)."""
    rng, lines = random.Random(5), []
    for i in range(n):
        start = rng.randrange(len(gene) - length)
        read = gene[start : start + length].encode()
        if strand == "reverse" or (strand == "both" and i % 2):
            read = reference.reverse_complement(read)
        lines += [f"@r{i}", read.decode(), "+", "I" * length]
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def dna_index(tmp_path: Path) -> Path:
    fasta = tmp_path / "genes.fna"
    fasta.write_text("".join(f">gene{i}\n{g}\n" for i, g in enumerate(GENES)))
    args = ["index", str(fasta), str(tmp_path / "idx"), "--targets", "--alphabet", "dna"]
    result = CliRunner().invoke(app, [*args, "--t-base", "1"])
    assert result.exit_code == 0, result.output
    return tmp_path / "idx"


def test_dna_index_records_its_molecule(dna_index: Path) -> None:
    meta = json.loads((dna_index / "meta.json").read_text())
    assert meta["molecule"] == "nucleotide"
    assert meta["params"]["alphabet"] == "dna"
    assert meta["params"]["k"] == DNA_K
    assert Index.load(dna_index).units.frame()["cluster_rep"].to_list() == [
        "gene0",
        "gene1",
        "gene2",
    ]


def test_reads_from_either_strand_hit_their_target(dna_index: Path, tmp_path: Path) -> None:
    reads = write_reads(tmp_path / "r.fastq", GENES[0], 200)
    result = CliRunner().invoke(
        app, ["query", str(dna_index), str(reads), "--out", str(tmp_path / "p.tsv"), "--aai"]
    )
    assert result.exit_code == 0, result.output
    got = pl.read_csv(tmp_path / "p.tsv", separator="\t")
    assert got["unit"].to_list() == [0]  # read hits only on gene 0
    assert got["containment"][0] > 0.95
    assert got["ani"][0] > 0.99
    assert not [c for c in got.columns if c.startswith("aai")]
    assert {"ani_naive", "ani_lo", "ani_hi"} <= set(got.columns)


def test_strand_does_not_matter(dna_index: Path, tmp_path: Path) -> None:
    index = Index.load(dna_index)
    forward = profile(index, write_reads(tmp_path / "f.fastq", GENES[1], 100, "forward"))
    reverse = profile(index, write_reads(tmp_path / "r.fastq", GENES[1], 100, "reverse"))
    assert forward.height == 1
    assert_frame_equal(forward, reverse)  # equal up to summation order


def test_nucleotide_and_protein_indexes_are_not_queried_jointly(
    dna_index: Path, tmp_path: Path
) -> None:
    proteins = write_members(tmp_path / "p.parquet", ["MKVLAAGIVLLAGHHKKERW" * 3])
    build_index(proteins, tmp_path / "aa", IndexParams(k=11, t_base=1.0))
    reads = write_reads(tmp_path / "r.fastq", GENES[0], 10)
    with pytest.raises(ValueError, match="nucleotide and protein"):
        profile(Index.load(dna_index), reads, extra=[Index.load(tmp_path / "aa")])
