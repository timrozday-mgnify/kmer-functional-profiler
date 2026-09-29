# mgnify-subset

Extracts whole MGnify90 clusters that have a member in a chosen biome from a MGnify
Proteins release, writes them as the index's input table, and optionally builds the
index for one or more nested subsets. It also runs the phase-6 cost study: per-cluster
statistics of the whole release and the index size they predict. Meant for HPC: the
release is about 1.3 TB of Parquet.

```text
MEMBERSHIP  mgy_clusters (biome filter, 1-in-N sample) + mgy_cluster_seqs -> membership/ (one
            directory per protein_id range)
EXTRACT     per protein_id range: mgy_protein_sequences, mgy_proteins_pfam  (--shards jobs)
MERGE       per cluster bucket -> members[.bucketB].parquet (protein_id, cluster_rep,
            full_length, sequence), pfam.parquet                          (--buckets jobs)
SUBSET      per sample N: clusters with cluster_rep % N == 0             (--sample list)
INDEX       kmer-functional-profiler index -> 1inN/index/               (--index false to skip)
STATS       per bucket: n_kmers, lengths per cluster; (hash, cluster) pairs at --stats_rate
COMBINE     -> clusters.parquet, groups.tsv, cost.tsv                     (--stats true)
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

### Cost study (phase 6)

Two runs, independent of each other.

Nested all-biome subsets, each built with the Python index; `trace.tsv` gives each
build's time and peak memory (`INDEX (1 in N)` rows), `1inN/index/meta.json` its sizes:

```bash
nextflow run workflows/mgnify-subset -profile slurm --biome root \
    --sample 100,1000,10000 --outdir cost-nested
```

Whole-release statistics. Every cluster lands in one of `--buckets` members tables
(about 4 GB each at 256), so each STATS job counts its clusters' distinct k-mers exactly,
with the build's own code (adapter masking, k, alphabet). COMBINE predicts the index for
each 1-in-N subset (`--cost_samples`) and each parameter set (`--cost_args`), so the
predictions for the nested subsets check against the builds above:

```bash
nextflow run workflows/mgnify-subset -profile slurm --biome root --sample 1 \
    --shards 256 --buckets 256 --pfam false --index false --stats true \
    --publish_mode link --outdir cost-full
```

This reads the whole sequence file once and writes it again as cluster buckets
(~1 TB under `cost-full/`, reused by the full build); use a local mirror if you have one.

| Output | Content |
| --- | --- |
| `clusters.parquet` | Per cluster: `cluster_rep`, `n_members`, `n_full_length`, `sum_len`, `max_len`, `n_kmers` |
| `groups.tsv` | Per sample: sampled k-mers (`hashes`) by the number of clusters holding them (`n_groups`) |
| `cost.tsv` | Per sample and parameter set: units, floored units, `below_floor` (floored units expecting fewer than `n_min` candidates), `t_max`, expected postings and tier/dense sizes in bytes |

Predicted sizes assume one key per posting (an upper bound) and `--sets` value sets
per posting (default 1, also an upper bound); the nested builds give the real ratio
(`tier2_sets / postings` in `meta.json`) to re-run COMBINE with.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--release` | EBI FTP `current_release` | Directory or https prefix holding the release Parquet files |
| `--biome` | `root:Host-associated:Human:Digestive system` | Substring of `cluster_members_biomes`; `root` for everything |
| `--sample` | `1000` | Keep clusters with `cluster_rep % sample == 0`; a list (`100,1000`) gives nested subsets, each a multiple of the smallest |
| `--shards` | `64` | Parallel `protein_id` ranges in EXTRACT |
| `--buckets` | `1` | Members tables, split by a hash of `cluster_rep`; `--index` needs 1 |
| `--max_protein_id` | `11200000000` | Upper end of the ranges; MEMBERSHIP fails if members lie above it |
| `--pfam` | `true` | Extract Pfam hits (needed by `--index`) |
| `--publish_mode` | `copy` | How members tables are published; `link` at full scale |
| `--index` | `true` | Build an index per sample, under `1inN/index/` |
| `--index_args` | `''` | Extra `kmer-functional-profiler index` options, e.g. `'--t-base 0.001 --n-min 8'` |
| `--stats` | `false` | Per-cluster statistics and predicted index sizes |
| `--stats_args` | `''` | `--k`, `--alphabet` for STATS and COMBINE |
| `--stats_rate` | `0.001` | Rate of the sampled (hash, cluster) pairs that measure promiscuity |
| `--cost_samples` | `1,10,100,1000,10000` | Subsets (`cluster_rep % N == 0`) COMBINE predicts |
| `--cost_args` | `--n-min 4 8 16 --t-cap 0.05 0.2 --t-dense 0 0.02 0.1` | Parameter grid for COMBINE (also `--t-base`, `--max-groups`, `--sets`) |
| `--venv` | repo `.venv` | Environment created by `setup.sh` |

Resources are set per label in `nextflow.config` (`process_medium` for DuckDB steps and
STATS, `process_high_memory` for INDEX and COMBINE); override them with `-c my.config`.
The index build holds the members table and the sampled k-mer tables in memory; the
nested run measures how its peak memory grows with the subset.
