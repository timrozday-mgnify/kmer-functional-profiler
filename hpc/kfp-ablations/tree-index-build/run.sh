#!/usr/bin/env bash
#SBATCH --job-name=kfp-tree-index-build
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Tree index, the three tiers on the 1-in-100 subset (docs/tree-index-plan.md, Progress log:
# tiers chosen from the per-protein survey), with workflows/tree-index/build.nf:
#   quick   t_base 0.01, n_min 8, tau 0, n_1 1
#   medium  t_base 0.02, n_min 16, tau 0.3, n_1 2
#   large   t_base 0.05, n_min 32, tau 0.3, n_1 2
# a leaf per distinct dense 8-mer set, the dense 8-mer topology. PARTITION (~130 GB, one job)
# splits ~53 M proteins into parts of 2 M; BUILD runs per part and tier (~26 GB each); MERGE
# writes one index per tier.
# Read:
#   results/index-<tier>/meta.json   postings by kind, bytes, positions, parts, build time
#   results/partition.json           components, the largest, part sizes
#   results/trace.tsv                time and peak RSS per task
# Re-submitting resumes.
source ../common.sh
[[ -f $SUBSET/members.parquet ]] || { echo "no $SUBSET/members.parquet" >&2; exit 1; }
nextflow run "$WF/tree-index/build.nf" -c ../site.config -resume --outdir results \
    --members "$(realpath "$SUBSET/members.parquet")"
