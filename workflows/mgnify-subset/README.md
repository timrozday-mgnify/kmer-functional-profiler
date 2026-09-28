# mgnify-subset

Extracts whole MGnify90 clusters that have a member in a chosen biome from a MGnify
Proteins release, writes them as the index's input table, and optionally builds the
index. Meant for HPC: the release is about 1.3 TB of Parquet.

```text
MEMBERSHIP  mgy_clusters (biome filter, 1-in-N sample) + mgy_cluster_seqs -> membership.parquet
EXTRACT     per protein_id range: mgy_protein_sequences, mgy_proteins_pfam  (--shards jobs)
MERGE       -> members.parquet (protein_id, cluster_rep, full_length, sequence), pfam.parquet
INDEX       kmer-functional-profiler index -> index/   (--index false to skip)
```

All members of a selected cluster are extracted, including members from other biomes,
so `p_in` is computed over the whole cluster.

## Setup

Needs Java 17+, Nextflow, [uv](https://docs.astral.sh/uv/) and a Rust toolchain.

```bash
bash workflows/setup.sh   # creates .venv in the repo root
```

## Run

Test profile (a 12-cluster fake release in `tests/data/mini_release`, a few seconds):

```bash
nextflow run workflows/mgnify-subset -profile test
```

Development subset (human gut, 1 in 1000 clusters) on Slurm:

```bash
nextflow run workflows/mgnify-subset -profile slurm --outdir gut-1in1000
```

With a local mirror of the release, point `--release` at it instead of the FTP site;
reading remotely works but moves the whole sequence file over HTTP.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--release` | EBI FTP `current_release` | Directory or https prefix holding the release Parquet files |
| `--biome` | `root:Host-associated:Human:Digestive system` | Substring of `cluster_members_biomes`; `root` for everything |
| `--sample` | `1000` | Keep clusters with `cluster_rep % sample == 0` |
| `--shards` | `64` | Parallel `protein_id` ranges in EXTRACT |
| `--max_protein_id` | `11200000000` | Upper end of the ranges; MERGE fails if members lie above it |
| `--index` | `true` | Build the index after extraction |
| `--index_args` | `''` | Extra `kmer-functional-profiler index` options, e.g. `'--t-base 0.001 --n-min 8'` |
| `--venv` | repo `.venv` | Environment created by `setup.sh` |

Resources are set per label in `nextflow.config` (`process_medium` for DuckDB steps,
`process_high_memory` for INDEX); override them with `-c my.config`. The index build
holds the members table and the sampled k-mer tables in memory; its peak memory has not
been measured on a real subset yet.
