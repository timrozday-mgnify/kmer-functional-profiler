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
uv run kmer-functional-profiler query idx/ r1.fq.gz r2.fq.gz --out profile.tsv
uv run kmer-functional-profiler import-sourmash KOs_sketched_scaled_1000.sig.zip ko_idx/ --ksize 11
nextflow run workflows/mgnify-subset -profile test   # MGnify subset + index; see its README
nextflow run workflows/fmh-benchmark -profile test   # benchmark vs fmh-funprofiler; see its README
```

`tests/python/test_compat.py` also checks parity with `sourmash prefetch` on
fmh-funprofiler's demo when its data is in `data/fmh-funprofiler/`:

```bash
mkdir -p data/fmh-funprofiler && cd data/fmh-funprofiler
curl -LO https://zenodo.org/records/10045253/files/KOs_sketched_scaled_1000.sig.zip
curl -LO https://raw.githubusercontent.com/KoslickiLab/fmh-funprofiler/main/demo/metagenome_example.fastq
```

```python
from kmer_functional_profiler import _core

for batch in _core.FastxHits("r1.fq.gz", "r2.fq.gz", k=11, max_hash=_core.max_hash(1 / 1000)):
    batch["read"], batch["mate"], batch["frame"], batch["hash"]  # numpy columns
```

Python prototype in `python/`, Rust kernels in `crates/core` (no PyO3),
bindings in `crates/py`.

Licence: GPL-3.0-or-later.
