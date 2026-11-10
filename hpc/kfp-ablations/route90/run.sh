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
# workflows/route90: MGnify90 units, k 9, keys by sqrt(members), 30-bit floor, confirmation on
# each unit's representative, every read pair recruited to the confirmed representatives (timed
# apart: not the index's cost). Three index sizes, each built and queried: 1 in 32 (standard),
# 1 in 16 and 1 in 8 (the larger options); RATES=32 (or any of them, comma-separated) for fewer.
# Read, per size <tag> (r32, r16, r8):
#   results/route90_recall_<tag>.tsv   per set x depth x AAI: target recall (the representative
#                                      within 0.9 of the best; Wilson interval), within 0.7 and
#                                      0.5, the unit's, confirmed share, read recall, precision
#   results/route90_targets_<tag>.tsv  present units, precision, strains lost before/after confirm
#   results/build_times_<tag>.json     postings, bytes per index file, seconds per shard
#   results/query_times_<tag>.json     seconds per stage (index, recruitment), units, peak RSS
#   results/trace.tsv                  CPU, wall and peak RSS per step
# Re-submitting resumes.
source ../common.sh
SLICE=../tree-index-slice/results
[[ -f $SLICE/members.parquet ]] || { echo "no $SLICE/members.parquet: run tree-index-slice first" >&2; exit 1; }
nextflow run "$WF/route90" -c ../site.config -resume --outdir results \
    --members "$(realpath "$SLICE/members.parquet")" --rates "${RATES:-32,16,8}"
