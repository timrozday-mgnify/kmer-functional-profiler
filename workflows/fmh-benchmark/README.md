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
                 overlaps one of its genes (the paper's rule, on CAMISIM's alignments)
PROFILE          kmer-functional-profiler query, every metagenome x every index
SCORE, SUMMARY   purity, completeness, completeness of the 25% least-covered true KOs,
                 base-weighted completeness per --min_hits value -> summary.tsv (mean, sd),
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
| `--min_hits` | `1,2,3,5` | Distinct k-mers for a KO to count as detected; each value is scored from the same profiles |
| `--indexes` | four configs (see `nextflow.config`) | `[name:, args:]` maps of `index` options |

Defaults compare, at k = 11: fmh-funprofiler's sketches (scaled 1000); our index at the
same base rate without and with the per-KO floor (`--n-min 8`); and our index at 10x
density (`--t-base 0.01`), the plan's "scaled = 100" baseline, without and with the floor. `index_*/meta.json` in the
output records each index's size; `trace.tsv` records each PROFILE task's runtime.
