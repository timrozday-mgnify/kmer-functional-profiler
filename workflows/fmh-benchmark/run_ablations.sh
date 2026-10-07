#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablations
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=kfp-ablations-%j.log
# Nextflow head job for one phase-7 ablation run (README, "Ablations"):
#   sbatch workflows/fmh-benchmark/run_ablations.sh index|reads|host|floor|tiers|standard|frames|aai|study OUTDIR [nextflow options]
# Run from the repo root, after workflows/setup.sh, with hpc.config next to this script. Each
# run has its own OUTDIR and work directory, so they can all run at once (fetch the shared
# databases first: README, "Ablations").
set -euo pipefail
here=workflows/fmh-benchmark
which=${1:?index, reads, host, floor, tiers, standard, frames, aai or study}
outdir=${2:?output directory}
shift 2
# module load nextflow singularity   # or apptainer; whatever your site provides
export NXF_OPTS='-Xms1g -Xmx4g'
nextflow run "$here" -profile slurm,singularity -c "$here/hpc.config" \
    -c "$here/ablations/$which.config" --outdir "$outdir" -w "$outdir/work" -resume "$@"
