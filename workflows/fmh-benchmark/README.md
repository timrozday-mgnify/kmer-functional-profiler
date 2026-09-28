# fmh-benchmark

KO detection by kmer-functional-profiler indexes against fmh-funprofiler's KO sketches on
simulated metagenomes, following the fmh-funprofiler paper (Bioinformatics 2024) but with
InSilicoSeq instead of CAMISIM. This is the phase-3 gate: completeness and purity at equal
density.

```text
FETCH            Zenodo 10055954 (CC-BY): KEGG genomes (9.1 GB zip), KEGG proteins (3.2 GB),
                 gene -> KO table, KO sketches. Downloaded once into --data_dir.
MEMBERS          proteins grouped by KO -> members.parquet
INDEX            kmer-functional-profiler index, one per --indexes entry
IMPORT_SKETCHES  the KO sketches as a sourmash-compatible index ("fmh_compat")
SAMPLE           --n_genomes random genomes per replicate (seed = replicate number)
SIMULATE         iss generate, lognormal abundances
TRUTH            reads mapped back with minimap2 (mappy); a KO is present if a read
                 overlaps one of its genes (the paper's rule, on CAMISIM's alignments); its
                 depth is aligned bases / gene length, summed over its genes
PROFILE          kmer-functional-profiler query, every metagenome x every index
SCORE, SUMMARY   purity, completeness, completeness of the 25% least-covered true KOs,
                 base-weighted completeness per count (kmers_hit; kmers_unique after gather;
                 kmers_wta, kmers_ufirst after winner-take-all, uniqueness-first)
                 and --min_hits value; abundance of the detected KOs against truth depth
                 (spearman_tp, l1) -> summary.tsv (mean, sd),
                 scores.tsv; profiles/ and truth/ keep the per-sample tables
```

The `fmh_compat` row reproduces fmh-funprofiler's KO calls: its hits equal the overlaps
`sourmash prefetch` reports (tested in `tests/python/test_compat.py`), and a KO counts as
detected at one shared hash, as with fmh-funprofiler's default `--threshold-bp 1000` at
scaled 1000. The paper's reference numbers (CAMISIM, wgsim, 64 genomes, 1 Gbp, k = 11):
purity 0.98, completeness 0.61 at 0% error; they are not directly comparable to runs
here, because the simulator differs.

## Setup

Needs Java 17+, Nextflow, [uv](https://docs.astral.sh/uv/) and a Rust toolchain.

```bash
bash workflows/setup.sh   # .venv with the package, InSilicoSeq and mappy
```

## Run

Test profile (tiny fake inputs in `tests/data/mini_fmh`, under a minute):

```bash
nextflow run workflows/fmh-benchmark -profile test
```

Full run on Slurm (about 13 GB download on first run, kept in `--data_dir`):

```bash
nextflow run workflows/fmh-benchmark -profile slurm --data_dir /path/to/fmh-benchmark-data
```

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--data_dir` | `fmh-benchmark-data` | Where the Zenodo inputs are stored (reused across runs) |
| `--sketches` | `KOs_sketched_scaled_1000.sig.zip` | KO sketches from the same Zenodo record |
| `--ksize` | `11` | Protein k of the sketches used |
| `--n_genomes` | `64` | Genomes per metagenome |
| `--replicates` | `10` | Metagenomes (seeds 1..N) |
| `--n_reads` | `6600000` | InSilicoSeq reads, both mates (~1 Gbp at 151 bp) |
| `--iss_model` | `novaseq` | InSilicoSeq error model: `hiseq`, `novaseq` or `miseq` |
| `--draws` | `100` | Posterior draws for 95% intervals on `coverage_zi` / `abundance_zi` and ambiguity groups (0: none) |
| `--min_hits` | `1,2` | Distinct k-mers for a KO to count as detected; each value is scored from the same profiles |
| `--indexes` | four configs (see `nextflow.config`) | `[name:, args:]` maps of `index` options |

Defaults compare, at k = 11: fmh-funprofiler's sketches (scaled 1000); our index at the
same base rate without and with the per-KO floor (`--n-min 8`); and our index at 10x
density (`--t-base 0.01`), the plan's "scaled = 100" baseline, alone and with a dense tier
at 1 in 10 (`--t-dense 0.1`: detection unchanged, EM abundances fitted on 10x more
k-mers). `index_*/meta.json` in the output records each index's size (`dense_bytes` for the
dense tier); `trace.tsv` records each PROFILE task's runtime.

With `-resume`, INDEX, IMPORT_SKETCHES and PROFILE rerun when the Python package changes,
and TRUTH and SCORE when `bench.py` does; Rust kernel changes need a fresh run.

Abundance is scored per row (`abundance` column) with an estimate that goes with its count:
`coverage` (hits per kept k-mer) with `kmers_hit`; `coverage_em` (EM over the units gather
keeps), `coverage_zi` (zero-inflated EM: coverage of the k-mers present) and `coverage_zib`
(zero-inflated with an empirical-Bayes prior on the present fraction) and `coverage_zip`
(zero-inflated with each k-mer's presence proportional to its in-KO frequency *p_in*) and
`abundance_zi` (`coverage_zi` times the gene copies present) with `kmers_unique`; and the hits each one-pass rule assigns per kept k-mer with its count:
`coverage_wta` with `kmers_wta` (each hit k-mer to the holding KO with the highest
containment, as sylph) and `coverage_ufirst` with `kmers_ufirst` (to the holding KO with the
highest Σ 1 / KOs-per-hit-k-mer, scaled by 1 / *t_g*). `spearman_tp` is the rank correlation with truth depth over true positives;
`l1` is the L1 distance between relative abundances over all true and detected KOs (0 is
exact, 2 is disjoint).

## Results (10 metagenomes, InSilicoSeq novaseq)

Means over seeds 1..10; sd of purity and completeness <= 0.015. "Low 25%" is completeness
on the least-covered quarter of true KOs. `kmers_hit` counts every hit k-mer for every KO
holding it; `kmers_unique` counts only the k-mers gather assigns to the KO.

| Index | Count | min_hits | Purity | Completeness | Low 25% |
| --- | --- | --- | --- | --- | --- |
| fmh_compat | kmers_hit | 1 | 0.975 | 0.688 | 0.295 |
| fmh_compat | kmers_unique | 1 | 0.991 | 0.672 | 0.268 |
| kfp_s1000 | kmers_hit | 1 | 0.985 | 0.677 | 0.285 |
| kfp_s1000 | kmers_unique | 1 | 0.996 | 0.663 | 0.263 |
| kfp_s1000_floor8 | kmers_hit | 1 | 0.976 | 0.719 | 0.404 |
| kfp_s100 | kmers_hit | 1 | 0.951 | 0.960 | 0.853 |
| kfp_s100 | kmers_hit | 2 | 0.976 | 0.911 | 0.694 |
| kfp_s100 | kmers_hit | 3 | 0.985 | 0.858 | 0.533 |
| **kfp_s100** | **kmers_unique** | **1** | **0.985** | **0.955** | **0.837** |
| kfp_s100 | kmers_unique | 2 | 0.994 | 0.903 | 0.665 |

`kfp_s100` scored on `kmers_unique` at `min_hits` 1 is the recommended setting: the purity
of `kmers_hit` at 3 with 10 points more completeness, and 27 points more than fmh_compat.
Gather removes 70% of the false positives (348 to 105 per sample; none left in all 10
samples, against 29 before) for 31 true KOs per sample, 90% of them in the lowest-coverage
quarter: single-k-mer hits taken by a more abundant relative (PTS and ABC transporter
paralogs, for example). The false positives it removes are mostly modular PKS/NRPS KOs
(pks2, pks8, tyrocidine and rapamycin synthases), whose shared domains carry identical
k-mers. Collagen VII (`K16628`) is the most persistent one left, in 5 of 10 samples.

At scaled 1000, `min_hits` above 1 costs most low-coverage KOs (completeness 0.49 at 2).
The floor helps at scaled 1000 but not at scaled 100, where nearly every KO already
samples more than 8 k-mers (`n_min 8` at `t_base 0.01` scored within 0.005 of `n_min 0`),
so that configuration was dropped. Gather adds no measurable query time (10-15 s per
sample either way).

### Abundance (same run design, `kmers_unique` >= 1, true positives)

| Index | Abundance | Spearman | L1 |
| --- | --- | --- | --- |
| fmh_compat | coverage | 0.08 | 1.41 |
| kfp_s100 | coverage | 0.32 | 1.35 |
| kfp_s100 | coverage_em | 0.35 | 1.25 |
| kfp_s100 | coverage_zi | 0.58 | 0.94 |
| kfp_s100_d10 | coverage_zi | 0.63 | 0.93 |
| kfp_s1000 | abundance_zi | 0.81 | 0.57 |
| kfp_s100 | **abundance_zi** | **0.96** | **0.22** |
| kfp_s100_d10 | abundance_zi | 0.99 | 0.14 |

Truth depth sums over each KO's gene copies in the sample, while `coverage_zi` is depth per
copy; `abundance_zi` scales it by the copies present (present k-mers / the k-mers an average
member holds) and matches truth closely. Recommended: `kfp_s100`, `kmers_unique` >= 1,
`abundance_zi`. For `fmh_compat` it reduces to plain EM, since imported sketches carry no
members. The dense tier costs 15x the table size, ~6x the query time and 36.5 GB to build.
