#!/usr/bin/env bash
#SBATCH --job-name=kfp-query-memory
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Query memory at full scale, as Slurm counts it (plan, phase 6, step 40): the release-tiers
# indexes queried in place (no scratch copy, no preload: as the pileup's QUERY runs), each
# cell at a ladder of memory limits. Needs the finished release-tiers run.
# Read results/query_cost.tsv, one row per (tier, pairs, limit):
#   end_exit            0, or 137 where the cgroup's OOM killer stopped the query at that limit
#   running             the stage a killed query died in
#   end_cgroup_peak     the job cgroup's high-water mark (what the limit is checked against)
#   {stage}_peak_cgroup and _peak_cgroup_{anon,file,mapped,shmem}: the charge and its parts
#                       by stage (file/mapped: index pages); {stage}_peak_anon, _peak_rss
# The smallest limit with end_exit 0 is what a query needs today; the 480 GB cells give the
# charge with nothing reclaimed. TIERS='slim standard large' sbatch run.sh adds the others.
source ../common.sh
TIERS=${TIERS:-standard}
for t in $TIERS; do
    [[ -f ../release-tiers/results-$t/index/meta.json ]] ||
        { echo "release-tiers has not finished (no results-$t/index)" >&2; exit 1; }
done
tiers=$(printf '%s,' $TIERS)
tiers=${tiers%,}
[[ $tiers == *,* ]] && tiers="{$tiers}"  # a one-name brace group matches nothing in Nextflow
nextflow run "$WF/mgnify-subset" -c ../site.config --query true --index false \
    --query_indexes "../release-tiers/results-$tiers/index" \
    --query_run ERR7738575,ERR7746321 --query_ladder 10000000,100000000 \
    --query_draws 0 --query_scratch '' --query_preload false \
    --query_memory '32 GB,64 GB,128 GB,192 GB,256 GB,480 GB' \
    --outdir results -resume
