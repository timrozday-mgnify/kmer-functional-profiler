#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-frames-query-cost
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# after frames-novaseq has finished: query cost of its three indexes in each frame mode
# (pooled gut sample, no draws), one outdir per mode: results-stopfree, results-edges20, ...
[[ -f ../frames-novaseq/results/summary.tsv ]] || { echo 'frames-novaseq has not finished' >&2; exit 1; }
for mode in stopfree edges:20 edges:30 all; do
    # a config file, not --query_args: Nextflow reads a value starting with '--' as an option
    printf 'params.query_args = "--frames %s"\n' "$mode" > "frames-${mode/:/}.config"
    nextflow run "$WF/mgnify-subset" -c ../site.config -c "frames-${mode/:/}.config" \
        --query true --index false \
        --query_indexes '../frames-novaseq/results/mgnify_index/*/index' \
        --query_run ERR7738575,ERR7746321 --query_ladder 4000000,40000000,100000000 \
        --query_draws 0 --outdir "results-${mode/:/}" -resume
done
