#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-standard-query-cost
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# after standard-novaseq has finished: query cost of its five indexes (pooled gut sample, no draws)
[[ -f ../standard-novaseq/results/summary.tsv ]] || { echo 'standard-novaseq has not finished' >&2; exit 1; }
nextflow run "$WF/mgnify-subset" -c ../site.config --query true --index false \
    --query_indexes '../standard-novaseq/results/mgnify_index/*/index' \
    --query_run ERR7738575,ERR7746321 --query_ladder 4000000,40000000,100000000 \
    --query_draws 0 --outdir results -resume
