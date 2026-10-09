#!/usr/bin/env bash
#SBATCH --job-name=kfp-aai-rescore
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# The floor runs' MGnify profiles rescored against nearest-member identity, aai re-estimated
# under the step-34 survival model and two comparisons (plan, phase 7, step 34). Needs the
# finished floor-novaseq and floor-1.65M runs; rebuilds no index and reruns no query.
# (The step34q check on the large tier, tbase05_all, came from abundance-tiers: plan, phase 7, step 38.)
source ../common.sh
for r in floor-novaseq floor-1.65M; do
    nextflow run "$WF/aai-rescore" -c ../site.config --run "../$r/results" \
        --members "$SUBSET/members.parquet" --genomes ../fmh-benchmark-data/genomes_extracted_from_kegg \
        --outdir "results/$r" -w "work/$r" -resume
done
