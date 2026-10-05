#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-tiers-query-cost
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# after tiers-novaseq has finished: query cost of its seven release-tier indexes (pooled gut sample, no draws)
[[ -f ../tiers-novaseq/results/summary.tsv ]] || { echo 'tiers-novaseq has not finished' >&2; exit 1; }
nextflow run "$WF/mgnify-subset" -c ../site.config --query true --index false \
    --query_indexes '../tiers-novaseq/results/mgnify_index/*/index' \
    --query_run ERR7738575,ERR7746321 --query_ladder 4000000,40000000,100000000 \
    --query_draws 0 --outdir results -resume
