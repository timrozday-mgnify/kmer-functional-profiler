#!/usr/bin/env bash
#SBATCH --job-name=fmh-benchmark
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=fmh-benchmark-%j.log
# Nextflow head job for the fmh benchmark with other tools (see README.md, "HPC").
#   sbatch workflows/fmh-benchmark/run_hpc.sh [extra nextflow options, e.g. --tools diamond]
# Run from the repo root, after workflows/setup.sh, with hpc.config next to this script.
set -euo pipefail
here=workflows/fmh-benchmark
# module load nextflow singularity   # or apptainer; whatever your site provides
export NXF_OPTS='-Xms1g -Xmx4g'
nextflow run "$here" -profile slurm,singularity -c "$here/hpc.config" -resume "$@"
