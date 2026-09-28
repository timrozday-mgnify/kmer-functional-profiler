"""Throughput of the read kernel (phase 1 gate: >= 1 M reads/min/thread).

Run with ``uv run pytest tests/bench``; not part of the default test run.
"""

import random
from pathlib import Path

import pytest
from pytest_benchmark.fixture import BenchmarkFixture

from kmer_functional_profiler import _core

N_READS = 200_000
GATE_READS_PER_MIN = 1_000_000


@pytest.fixture(scope="module")
def reads(tmp_path_factory: pytest.TempPathFactory) -> Path:
    rng = random.Random(0)
    path = tmp_path_factory.mktemp("bench") / "reads.fastq"
    with path.open("w") as f:
        for i in range(N_READS):
            f.write(f"@r{i}\n{''.join(rng.choices('ACGT', k=150))}\n+\n{'I' * 150}\n")
    return path


@pytest.mark.parametrize("frames", ["stopfree", "all"])
def test_fastx_throughput(benchmark: BenchmarkFixture, reads: Path, frames: str) -> None:
    def run() -> int:
        return sum(
            len(b["hash"])
            for b in _core.FastxHits(reads, k=10, frames=frames, max_hash=_core.max_hash(0.01))
        )

    benchmark.pedantic(run, rounds=3)
    reads_per_min = N_READS / benchmark.stats.stats.mean * 60
    benchmark.extra_info["reads_per_min"] = reads_per_min
    assert reads_per_min >= GATE_READS_PER_MIN
