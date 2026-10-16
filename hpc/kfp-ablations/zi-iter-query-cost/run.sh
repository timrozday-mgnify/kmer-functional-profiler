#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-zi-iter-query-cost
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# The zero-inflated fit's step cap (plan, phase 7, step 43): promiscuity-novaseq's tbase05_all_g16
# (the large tier at the new max_groups) on the pooled gut sample at 40, 100 and 200 M pairs, with
# --zi-max-iter 1,000 (step 42's cap), 3,000 and 10,000 (the new default), one outdir per cap:
# results-zi1000, ... Read zi_em_iterations / _unconverged_units / _max_change and fit_zi_wall_s
# per cap (cost and convergence only: the query run publishes stats, not profiles; g16-novaseq
# measures the cap's effect on accuracy).
[[ -f ../promiscuity-novaseq/results/summary.tsv ]] || { echo 'promiscuity-novaseq has not finished' >&2; exit 1; }
for cap in 1000 3000 10000; do
    # a config file, not --query_args: Nextflow reads a value starting with '--' as an option
    printf 'params.query_args = "--zi-max-iter %s"\n' "$cap" > "zi${cap}.config"
    nextflow run "$WF/mgnify-subset" -c ../site.config -c "zi${cap}.config" \
        --query true --index false \
        --query_indexes '../promiscuity-novaseq/results/mgnify_index/tbase05_all_g16/index' \
        --query_run ERR7738575,ERR7746321 --query_ladder 40000000,100000000,200000000 \
        --query_draws 0 --outdir "results-zi${cap}" -resume
done
