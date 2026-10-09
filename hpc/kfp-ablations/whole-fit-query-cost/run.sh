#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-whole-fit-query-cost
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# Whole fits instead of block-wise (plan, phase 7, step 41): promiscuity-novaseq's two controls
# at 100 and 200 M pairs with --max-fit-pairs 100 M, so no component is fitted block-wise.
# Compare with promiscuity-query-cost's same cells (block-wise above 10 M pairs): fit_em and
# fit_zi wall, fit_*_peak_anon, em_* / zi_em_* (block rounds, unconverged units), and the profiles.
[[ -f ../promiscuity-novaseq/results/summary.tsv ]] || { echo 'promiscuity-novaseq has not finished' >&2; exit 1; }
nextflow run "$WF/mgnify-subset" -c ../site.config --query true --index false \
    --query_indexes '../promiscuity-novaseq/results/mgnify_index/tbase{05,10}_all/index' \
    --query_run ERR7738575,ERR7746321 --query_ladder 100000000,200000000 \
    --query_args '--max-fit-pairs 100000000' --query_draws 0 --outdir results -resume
