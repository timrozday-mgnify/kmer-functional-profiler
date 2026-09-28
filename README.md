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
```

Python prototype in `python/`, Rust kernels in `crates/core` (no PyO3),
bindings in `crates/py`.

Licence: GPL-3.0-or-later.
