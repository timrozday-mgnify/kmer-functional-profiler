#!/usr/bin/env bash
#SBATCH --job-name=kfp-route90
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# The tree index's route on its own on the compact index (docs/tree-index-plan.md, Progress log:
# the compact route, 2026-11-10), on the tree-index slice (4.2 M proteins), with
# workflows/route90: MGnify90 units, k 9, 1 in 32 by sqrt(members), 30-bit floor, confirmation
# on each unit's representative, every read pair recruited to the confirmed representatives.
# Read:
#   results/route90_recall.tsv    per set x depth x AAI: target recall (the representative within
#                                 0.9 of the best; Wilson interval), within 0.7 and 0.5, the
#                                 unit's, confirmed share, read recall and precision
#   results/route90_targets.tsv   present units, precision, strains lost before / after confirm
#   results/build_times.json      postings, bytes per index file, seconds per shard
#   results/query_times.json      seconds per stage, units seeded and confirmed, peak RSS
#   results/trace.tsv             CPU, wall and peak RSS per step
# Re-submitting resumes.
source ../common.sh
SLICE=../tree-index-slice/results
[[ -f $SLICE/members.parquet ]] || { echo "no $SLICE/members.parquet: run tree-index-slice first" >&2; exit 1; }
nextflow run "$WF/route90" -c ../site.config -resume --outdir results \
    --members "$(realpath "$SLICE/members.parquet")"
