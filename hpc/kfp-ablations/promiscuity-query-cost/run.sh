#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-promiscuity-query-cost
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# after promiscuity-novaseq: its seven indexes queried on the pooled gut sample at 40, 100 and
# 200 M pairs. Read detected_largest_component_*, fit_em/fit_zi wall, em_* and zi_em_* (the
# zero-inflated fit's convergence, new in step 40) and holders_cut_* (what each cut predicts
# from the 64 control). 200 M pairs is where tbase05_all may first go block-wise.
[[ -f ../promiscuity-novaseq/results/summary.tsv ]] || { echo 'promiscuity-novaseq has not finished' >&2; exit 1; }
nextflow run "$WF/mgnify-subset" -c ../site.config --query true --index false \
    --query_indexes '../promiscuity-novaseq/results/mgnify_index/*/index' \
    --query_run ERR7738575,ERR7746321 --query_ladder 40000000,100000000,200000000 \
    --query_draws 0 --outdir results -resume
