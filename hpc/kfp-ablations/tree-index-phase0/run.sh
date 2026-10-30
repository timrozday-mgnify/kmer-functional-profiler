#!/usr/bin/env bash
#SBATCH --job-name=kfp-tree-index-phase0
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Tree index, phase 0 (docs/tree-index-plan.md), with workflows/tree-index: the release's
# 256 members buckets (release-tiers' standard run) read twice. SAMPLE takes 1 in 5,000
# clusters (~330 k, ~1.2 M proteins). STREAM counts set digests, signatures and k-mer slices
# over the whole release and finds the representatives sharing >= 4 stored k-mers with the
# sample. MEASURE builds the trees on the sample (256 GB, up to 48 h), and TIER1 builds today's
# standard tier on the same proteins.
# Read:
#   results/gate.tsv                  the gate: postings <= 1.5x tier 1, leaves <= 2x units
#   results/measure-core/measure.json everything measured on the sample; TSVs beside it
#   results/measure-all/              the same with the cross-cluster representatives
#   results/stream.json               whole-release counts (equal-set collapse, signature
#                                     groups, k-mer holders and distinct/total per slice)
source ../common.sh
STD=../release-tiers/results-standard
compgen -G "$STD/members.bucket*.parquet" >/dev/null || { echo "no $STD/members.bucket*.parquet: run release-tiers first" >&2; exit 1; }
nextflow run "$WF/tree-index" -c ../site.config -resume --outdir results \
    --store "$(realpath "$STD")/members.bucket*.parquet" --sample 5000
