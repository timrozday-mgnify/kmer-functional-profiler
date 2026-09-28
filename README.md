# kmer-functional-profiler

Protein k-mer functional profiling of metagenomes: frame-filtered translated
reads, FracMinHash sampling with nested per-family thresholds, and an EM
abundance model over MGnify Proteins clusters. Early development; see
[the plan](docs/kmer-functional-profiler-plan.md).

## Development

```bash
uv sync                  # builds the Rust extension (release profile)
uv run pytest
cargo test --workspace
prek install --hook-type pre-commit --hook-type pre-push
uv run pytest tests/bench                        # throughput (phase 1 gate: >= 1 M reads/min/thread)
uv run python scripts/make_fixtures.py           # regenerate tests/data
uv run python scripts/fetch_mgnify_sample.py     # small MGnify sample into data/ (gitignored)
uv run kmer-functional-profiler index members.parquet idx/ --pfam pfam.parquet
nextflow run workflows/mgnify-subset -profile test   # MGnify subset + index; see its README
```

```python
from kmer_functional_profiler import _core

for batch in _core.FastxHits("r1.fq.gz", "r2.fq.gz", k=11, max_hash=_core.max_hash(1 / 1000)):
    batch["read"], batch["mate"], batch["frame"], batch["hash"]  # numpy columns
```

Python prototype in `python/`, Rust kernels in `crates/core` (no PyO3),
bindings in `crates/py`.

Licence: GPL-3.0-or-later.
