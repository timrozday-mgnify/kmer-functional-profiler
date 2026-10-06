#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-em-query-cost
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# EM convergence (plan, phase 7, step 29): tiers-novaseq's tbase05_all and tbase10_all indexes
# queried again with the new EM, in a fresh run directory so no QUERY is served from the cache.
# tbase10_all at 100 M pairs has the one component above MAX_FIT_PAIRS (block-wise; 1,777 s
# in the EM before); compare fit_em_wall_s, em_* and the profiles with tiers-query-cost's.
[[ -d ../tiers-novaseq/results/mgnify_index/tbase10_all/index ]] || { echo 'tiers-novaseq indexes missing' >&2; exit 1; }
nextflow run "$WF/mgnify-subset" -c ../site.config --query true --index false \
    --query_indexes '../tiers-novaseq/results/mgnify_index/tbase{05,10}_all/index' \
    --query_run ERR7738575,ERR7746321 --query_ladder 40000000,100000000 \
    --query_draws 0 --outdir results -resume
