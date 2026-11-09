#!/usr/bin/env bash
#SBATCH --job-name=kfp-route90
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# The tree index on its own, no anchoring (docs/tree-index-plan.md, Progress log: the route on
# its own, 2026-11-09), on the tree-index slice (4.2 M proteins), with workflows/route90:
# MGnify90 units, k 9, 1 in 8 by sqrt(members), 30-bit floor, confirmation on each unit's
# representative, pseudo-alignment to the confirmed units' leaves, seed alignment of the top hits.
# Read:
#   results/route90_recall.tsv    per set x depth x AAI: target recall (Wilson interval over
#                                 strains), confirmed share, read recall and precision
#   results/route90_targets.tsv   present targets, precision, strains lost before / after confirm
#   results/build_times.json      postings, bytes per index file, build time per part
#   results/query_times.json      seconds per stage and alignment lines per read set
#   results/trace.tsv             CPU, wall and peak RSS per step
# Re-submitting resumes.
source ../common.sh
SLICE=../tree-index-slice/results
[[ -f $SLICE/members.parquet ]] || { echo "no $SLICE/members.parquet: run tree-index-slice first" >&2; exit 1; }
nextflow run "$WF/route90" -c ../site.config -resume --outdir results \
    --members "$(realpath "$SLICE/members.parquet")"
