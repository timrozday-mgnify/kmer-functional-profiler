# aai-rescore

Rescores a finished fmh-benchmark run's MGnify unit profiles against **nearest-member
identity** (plan, phase 7, step 34), with `aai` re-estimated under survival models. Nothing is
rebuilt or re-queried: `aai-score --model` recomputes `aai` from each profile's saved
zero-inflated fit. That makes this the cheap way to score runs that finished before the
benchmark gained its nearest-member pass. New fmh-benchmark runs do the same themselves
(`--aai_models`).

```text
MEMBERS          members of every cluster in the run's mgnify/<annotation>_gene_units.parquet
PROTEINS         the benchmark genomes' proteins, in --chunks FASTA files
MGNIFY_DB, MGNIFY_ANNOTATE
                 DIAMOND blastp of those proteins against the members
NEAREST          per (gene, cluster) the nearest member -> results/gene_units.parquet
SCORE, SUMMARY   aai-score of every units/seed<sid>_<index>[~arm].tsv, as profiled and under
                 each --models entry (arm '<arm>+<name>') -> aai_summary.tsv, aai_scores.tsv
```

## Run

Needs the repo's venv (`bash workflows/setup.sh`), Nextflow, and Singularity or Docker
(DIAMOND). `hpc/kfp-ablations/aai-rescore/run.sh` runs it for the `floor` runs:

```bash
nextflow run workflows/aai-rescore -profile slurm,singularity --run <fmh outdir> \
    --members <subset>/members.parquet --genomes <data_dir>/genomes_extracted_from_kegg
```

`--models` defaults to the step-34 fit and two comparisons
(`workflows/fmh-benchmark/ablations/aai_models/`): `step34`, `step34_nounion` (union term
0) and `step34_ends1` (no end loss). `-stub -profile test` checks the wiring.

## Read

In `aai_summary.tsv`, `aai_truth` is `nearest`. Step 25's decision rule is bias against the
union truth within ±0.02 at 70–95% for units with `aai_kmers` ≥ 5, and interval coverage
≥ 0.9 per `aai_kmers` bin. Apply it to the `+step34` arms. `aai_cover_members<lo>` gives
coverage by the unit's members, which tests the union term's variance shortcut.
