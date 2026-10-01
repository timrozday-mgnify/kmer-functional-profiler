# mgnify-subset

Extracts whole MGnify90 clusters that have a member in a chosen biome from a MGnify
Proteins release, writes them as the index's input table, and optionally builds the
index for one or more nested subsets. It also runs the phase-6 cost study (per-cluster
statistics of the whole release and the index size they predict) and the partitioned
index build over cluster buckets. Meant for HPC: the
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
CANDIDATES  per bucket: unit table, candidate hashes                       (--build true)
BLOOM       all candidate hashes -> Bloom filter, t_max, hash quantiles
PRESENCE    per bucket: (hash, cluster) at t_max passing the filter, by hash range
GROUPS      per hash range: n_groups; candidate rows by bucket              (--ranges jobs)
POSTINGS    per bucket: score, promiscuity cut, floor, per-unit columns, Pfam
UNITS       all units numbered by cluster_rep; tier-2 layout and pack ranges
PACK_RANGE  per pack range: keys and value sets                            (--ranges jobs)
DEDUP       per set-hash range: distinct value sets                        (--ranges jobs)
CONCAT      -> index/
FETCH_READS paired ENA run(s) (--query_run)                                 (--query true)
POOL        several runs concatenated into one sample                  (two or more runs)
LADDER      nested read subsets of it                                    (--query_ladder)
QUERY       per subset x index: kmer-functional-profiler query --stats
QUERY_COST  -> query_cost.tsv
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

### Partitioned build (phase 6)

The single-process INDEX holds every member and k-mer table at once, which does not fit
the whole release on one node. `--build true` splits the same build
(`kmer_functional_profiler.partition`) over the MERGE buckets: everything a unit needs is
local to its bucket except `n_groups` (units holding a candidate k-mer), which a
hash-range reduce counts exactly. A blocked Bloom filter of all candidate hashes
(Rust, one cache line per key; `--bloom_bits` per hash, ~1% false hits at 10) keeps each bucket's rows for that reduce
to about 1.1 per candidate; GROUPS drops the false hits, so the index equals INDEX's on
the same members (checked by `tests/python/test_partition.py` and the `test_build`
profile). Hash ranges are cut at quantiles of the candidate hashes, so they hold about
equal numbers of candidates. No dense tier yet (`--t-dense` must be 0).

```bash
nextflow run workflows/mgnify-subset -profile test,test_build   # two buckets, a few seconds
nextflow run workflows/mgnify-subset -profile slurm --biome root --sample 1 \
    --shards 256 --buckets 256 --ranges 256 --index false --build true \
    --publish_mode link --outdir full-build
```

Tier 2 is packed in hash ranges too: UNITS fixes the table layout from the total
postings and cuts pack ranges at key boundaries, each PACK_RANGE packs its range's
postings, DEDUP removes duplicate value sets across ranges (split by a hash of the
sets' content), and CONCAT joins the parts and sets; sets are numbered in hash order, so
the result is the same as packing at once. No step holds every posting, and CONCAT
holds only the distinct sets.

At 1 in 1000 (8 buckets, 8 ranges, run serially on a laptop) the stages took 11 s
(CANDIDATES), 0.2 s (BLOOM), 10 s (PRESENCE), 1 s (GROUPS), 7 s (POSTINGS), 0.1 s
(UNITS), 1.2 s (PACK_RANGE), 0.7 s (DEDUP) and 0.2 s (CONCAT), with presence files of
~11.5 bytes per
row. Scaled to the whole release: a ~14 GB Bloom filter (about a minute to fill;
memory-mapped by every PRESENCE job, so jobs on one node share it), ~140 GB
of presence files, a ~34 GB unit table in UNITS, and CONCAT holding tier 2 (~31 GB).

### Query cost (phase 6, Q1)

Profiles nested read subsets of one real metagenome against each index given, with
the current query plus stage timers, to find which stages cost time and memory at scale.
The default run is ERR7738575 (human gut, NovaSeq, 111.5 M pairs = 223 M reads, 31 Gbp),
which LADDER subsets to 0.01, 0.1, 0.4, 1.2, 4, 12, 40 and 100 M pairs (100 M pairs is
200 M reads, about the deepest metagenome to expect); each subset contains the smaller
ones. Several comma-separated runs are fetched and pooled into one sample in the order
given, for depths and diversity beyond one run: e.g. `--query_run ERR7738575,ERR7746321`
(two people's gut, 217.6 M pairs) with `--query_ladder ...,100000000,200000000`. Step 31
of the plan used ERR7746321 (Hadza gut, 106 M pairs) up to 40 M pairs.

```bash
nextflow run workflows/mgnify-subset -profile test               # builds results-test/1in*/index
nextflow run workflows/mgnify-subset -profile test,test_query    # a few seconds
nextflow run workflows/mgnify-subset -profile slurm --query true --index false \
    --query_indexes 'cost-nested/1in*/index' --outdir query-cost
```

Compute nodes without internet: run FETCH_READS on the head node
(`process.withName: 'FETCH_READS' { executor = 'local' }` in a `-c` config), or download
the run and pass `--query_reads R1,R2`. Each QUERY writes `query/{index}.{pairs}.{draws}.json`
(`kmer-functional-profiler query --stats`); `query_cost.tsv` has one row per query:
the counts (reads, sampled and distinct sampled k-mers, hit k-mers, hit rows, (unit, hash)
pairs, hit and detected units, component sizes; `fit_batches` and
`fit_largest_batch_pairs` for the EM's batches of components, `em_iterations` and
`em_unconverged_units` for its convergence) and `{stage}_wall_s`, `_cpu_s` and
`_peak_rss` (bytes, the process's peak at the stage's end, including resident pages of
the memory-mapped index) and `_peak_anon` (the same peak for anonymous memory only, Linux
`RssAnon` sampled every 50 ms) per stage. Queries that exceed
`--query_memory` or 24 h fail without stopping the rest; they are FAILED in `trace.tsv`.

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
| `--index_args` | `''` | Extra `kmer-functional-profiler index` options, e.g. `'--t-base 0.001 --n-min 8'`; also used by `--build` |
| `--build` | `false` | Partitioned index build over the buckets, under `index/` (one `--sample`, `--pfam true`) |
| `--ranges` | `64` | Hash ranges of the partitioned build's `n_groups` reduce and tier-2 packing |
| `--bloom_bits` | `10` | Bloom filter bits per candidate hash |
| `--stats` | `false` | Per-cluster statistics and predicted index sizes |
| `--stats_args` | `''` | `--k`, `--alphabet` for STATS and COMBINE |
| `--stats_rate` | `0.001` | Rate of the sampled (hash, cluster) pairs that measure promiscuity |
| `--cost_samples` | `1,10,100,1000,10000` | Subsets (`cluster_rep % N == 0`) COMBINE predicts |
| `--cost_args` | `--n-min 4 8 16 --t-cap 0.05 0.2 --t-dense 0 0.02 0.1` | Parameter grid for COMBINE (also `--t-base`, `--max-groups`, `--sets`) |
| `--query` | `false` | Query cost study (skips extraction unless `--index`, `--stats` or `--build`) |
| `--query_indexes` | `''` | Index directories (glob), each named by its parent directory |
| `--query_run` | `ERR7738575` | Paired ENA run to download; comma-separated runs are pooled into one sample |
| `--query_reads` | `''` | `R1,R2` local files instead of `--query_run` |
| `--query_ladder` | `10000,100000,400000,1200000,4000000,12000000,40000000,100000000` | Nested subset sizes, in pairs (each at most the run's pairs) |
| `--query_seed` | `1` | Seed of the ladder's shuffle |
| `--query_draws` | `100` | Posterior draws, on the `--query_draws_on` index only (0 elsewhere) |
| `--query_draws_on` | `1in100` | Index that also gets the posterior |
| `--query_memory` | `128 GB` | Memory per query |
| `--query_scratch` | `${TMPDIR:-/tmp}` | Node-local directory to copy the index's query files to before each query (removed after; `false` queries it in place). If the copy fails, the query uses the index in place |
| `--query_in_memory` | `false` | Read the tier-2 (and dense) arrays into memory instead of memory-mapping them: no page faults in lookup, but their full size (~31 GB at full scale) adds to peak RSS, so raise `--query_memory` |
| `--query_low_memory` | `false` | Trade time for memory where results are unchanged: with draws, the reads are hashed a second time for the detected units' k-mers instead of keeping per-read rows of every hit unit (~2 GB at 40 M pairs, one more read pass). For memory-capped runs; on HPC raising `--query_memory` is cheaper |
| `--query_preload` | `true` | Read tier 2 into the page cache before each query, so lookups are not random reads from a network filesystem (time logged to `.command.err`) |
| `--venv` | repo `.venv` | Environment created by `setup.sh` |

Resources are set per label in `nextflow.config` (`process_medium` for DuckDB steps and
STATS and the partitioned build's per-bucket and per-range steps, `process_high_memory`
for INDEX, COMBINE, UNITS and CONCAT); override them with `-c my.config`.
The index build holds the members table and the sampled k-mer tables in memory; the
nested run measures how its peak memory grows with the subset.
