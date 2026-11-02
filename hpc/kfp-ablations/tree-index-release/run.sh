#!/usr/bin/env bash
#SBATCH --job-name=kfp-tree-index-release
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Tree index on the whole MGnify Proteins release (docs/tree-index-plan.md, phase 5), with
# workflows/tree-index/build.nf from release-tiers' 256 members buckets (~5.7 x 10^9 proteins).
#   PARTITION  per bucket (--partition_per_file; each bucket holds whole MGnify90 clusters):
#              k 11 storage components packed into parts of 1 M proteins (~22 M proteins and
#              ~55 GB per bucket, from the 1-in-100 subset's ~130 GB for 53 M). Edges between
#              clusters in different buckets are lost: their k-mers are posted once per bucket.
#   BUILD      per part and tier (~5,700 parts per tier, ~30 GB each; 96 GB given)
#   MERGE      per tier, out of core by hash range (merge.py): -> results/index-<tier>/
#              (standard ~195-215 GB, large ~360-415 GB, slim ~100-130 GB; estimates)
#   STORE      per part, -> results/store-<tier>/ (flat, ~1x the residues: ~1.6 TB per tier)
# No benchmark here (--bench false): tree-index-reads and the phase-4 benchmark use the result.
# TIERS='standard:--mode standard' by default; add slim or large once tree-index-build has
# chosen. STORE_DIR moves the store into results/ (no second copy of ~1.6 TB; nothing reads it after).
# Read:
#   results/b*-partition.json       per bucket: proteins, components, largest, parts
#   results/index-<tier>/meta.json  postings by kind, bytes, positions, parts, build and merge
#                                   time
#   results/trace.tsv               time and peak RSS per task
# Re-submitting resumes.
source ../common.sh
STD=../release-tiers/results-standard
compgen -G "$STD/members.bucket*.parquet" >/dev/null || { echo "no $STD/members.bucket*.parquet: run release-tiers first" >&2; exit 1; }
cat > release.config <<'CFG'
process {
    withName: '.*PARTITION' { memory = { 96.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: '.*BUILD' { memory = { 96.GB * task.attempt } }
    withName: '.*MERGE' { cpus = 4; memory = { 128.GB * task.attempt }; time = { 72.h * task.attempt } }
    withName: '.*STORE_DIR' { publishDir = [path: 'results', mode: 'move'] }
}
CFG
nextflow run "$WF/tree-index/build.nf" -c ../site.config -c release.config -resume --outdir results \
    --members "$(realpath "$STD")/members.bucket*.parquet" --partition_per_file true --bench false --part_proteins 1000000 \
    --tiers "${TIERS:-standard:--mode standard}"
