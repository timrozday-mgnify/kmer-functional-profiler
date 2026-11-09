# pfam-reannotate

Searches MGnify proteins with `hmmsearch` against Pfam-A and keeps matches below the
gathering thresholds (GA). The release's own Pfam table (`mgy_proteins_pfam`) is Pfam 38.2
at `--cut_ga`, so it holds no weak hits. The default HMM file is the same Pfam version.

```text
PFAM_HMM   Pfam-A -> Pfam-A.hmm, ga.tsv (each model's GA thresholds)
PLAN       ids of the targets -> ranges.tsv, protein_id ranges of ~--chunk_proteins
FASTA      per range: the sequences
HMMSEARCH  per range: --domtblout and alignments, gzipped
PARSE      per range -> hits/chunkN.parquet
SUMMARY    -> summary.tsv, one row per arm
COMPARE    -> compare.tsv, real against reversed                      (--reversed true)
```

## Setup

Needs Java 17+, Nextflow and [uv](https://docs.astral.sh/uv/). `HMMSEARCH` runs HMMER 3.4 from
a container on HPC, or `hmmsearch` from the `PATH` without one.

```bash
bash workflows/setup.sh   # creates .venv in the repo root
```

## Run

Test profile (a few seconds, six mini HMMs):

```bash
nextflow run workflows/pfam-reannotate -profile test
```

On HPC, 1 in 100 representatives first (`hpc/kfp-ablations/pfam-reannotate/run.sh` does this
with the site config):

```bash
nextflow run workflows/pfam-reannotate -profile slurm,singularity --sample 100 --outdir pfam-1in100
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--targets` | `reps` | `reps`: MGnify90 representatives; `all`: every protein |
| `--sample` | `1000` | 1 in N targets by `protein_id`. With `all`, random proteins, not whole clusters |
| `--reversed` | `false` | Also search every sequence reversed, the null arm; twice the compute |
| `--pp` | `false` | Store each alignment's posterior probabilities |
| `--release` | FTP site | A local mirror avoids moving sequences over HTTP |
| `--hmm` | Pfam 38.2 | Local or remote `Pfam-A.hmm[.gz]` |
| `--chunk_proteins` | `500000` | Proteins per `hmmsearch` job |
| `--hmmsearch_args` | see config | Reporting thresholds |

## Expected time

Measured on a laptop on 10,000 MGnify representatives (168 aa on average): 225, 180 and 173 s
at `--cpu` 1, 2 and 4. CPU time stays near 345 core-s, so jobs use `--cpu 2`. Reporting
thresholds and writing alignments do not change the time. Loading Pfam-A costs ~10 s per job. That gives
1.07×10⁻⁴ job-seconds and 2.14×10⁻⁴ allocated core-seconds per residue:

| Targets | Residues | Jobs | Core-hours | Wall at 75 jobs at once |
| --- | ---: | ---: | ---: | ---: |
| 1 in 100 representatives, with `--reversed` | 2 × 2.5×10⁹ | 66 | ~300 | ~3 h |
| All representatives | 2.5×10¹¹ | ~3,300 | ~14,900 | ~4 days |
| All proteins | 1.07×10¹² | ~11,500 | ~63,600 | ~18 days |

These are at laptop core speed, and HPC cores may be slower. The 1-in-100 run's `trace.tsv`
gives the real rate to rescale them. Each job takes about 2.2 h at laptop speed. Wall time
scales with `executor.queueSize` (75 in the site config).

## Output

`hits/{real,reversed}.chunkN.parquet`, one row per domain, sorted by `protein_id`:
`protein_id`, `pfam_accession`, `seq_evalue`, `seq_score`, `c_evalue`, `i_evalue`, `score`,
`hmm_from`, `hmm_to`, `env_from`, `env_to`, `domain`, `ali_from`, `ali_to`, `cigar`, `pp` (with
`--pp`), `above_ga`. Compare by bit score: E-values here use Z = 10⁶ (E = 10⁶ × P).

Alignments are stored as operations, not residues:

- `cigar`: `M` a model position aligned to a residue, `I` a residue with no model position, `D`
  a model position with no residue.
- The aligned protein rebuilds from the sequence, `ali_from` and `cigar`: `M` and `I` take the
  next residues, and `D` is a gap. In the `reversed` arm the sequence is the reversed one, and
  `ali_from` counts along it.
- The aligned model rebuilds from the HMM's consensus, `hmm_from` and `cigar`.
- `pp`: HMMER's posterior probability line, one character per alignment column (`0`-`9`, `*`,
  and `.` at a gap).

Measured on 5,000 proteins: 148 bytes per protein with `pp` and 101 without, so ~5 GB for the
1-in-100 run with `pp` and its null arm (if the reversed arm were as hit-rich, an upper bound), and
~0.17 TB for all representatives without `pp`.

`work/` holds each chunk's FASTA and gzipped HMMER output until deleted, about 1.7 KB per
protein on the sample: ~55 GB for the 1-in-100 run with its null arm, up to ~2.8 TB for all
representatives. Delete it after the run.

## Reading the results

- **Chance hits.** `-E 10` at Z = 10⁶ over 30,134 models expects about 0.3 chance sequence
  hits per protein, so much of the weak set may be chance. `--reversed true` measures it:
  `compare.tsv` has domains per 5-bit bin of `score` from `score_from`, each side of GA, in the
  `real` and `reversed` arms, and `reversed / real` estimates the chance share of a bin.
  Reversal keeps length and composition and so is a mild null: repeats and low-complexity
  sequence can still hit, while real homologues of reversed sequences would inflate it.
- **`above_ga` at the boundary.** domtblout rounds scores to 0.1 bit, so a domain just under GA
  can read as above it: 8 of 3,261 on 5,000 proteins.
- **Filters.** hmmsearch's default filters stay on, so hits they reject are never reported;
  `--max` turns them off at a far higher cost (`--hmmsearch_args`).
