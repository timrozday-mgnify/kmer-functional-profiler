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
NEAREST          per (gene, cluster) the nearest member, among hits at query and subject
                 coverage >= --nearest_min_cov (0.8, as aai-model's pairs) -> results/gene_units.parquet
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

`--models` defaults to the step-34 fit, its refit on the query's measure and two comparisons
(`workflows/fmh-benchmark/ablations/aai_models/`): `step34`, `step34q` (fitted on
`shared` / `pin_sum` with the union scatter), `step34_nounion` (union term 0) and
`step34_ends1` (no end loss). `-stub -profile test` checks the wiring.

## Read

In `aai_summary.tsv`, `aai_truth` is `nearest`. Step 25's decision rule is read on units one
present gene hits (`aai_*_single_*`; plan, step 34): bias within ±0.02 at 70–95% and
interval coverage ≥ 0.9 per true-identity bin. By `aai_kmers` bin, coverage near 0.9 is
what a correct model gives (binning on hits conditions on the outcome). The union metrics
score strain mixes against an independent-union truth that overstates them. Read it on the
large tier where the run has it. `aai_cover_members<lo>` gives coverage by the unit's
members.
