#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-abundance-tiers
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# Abundance under divergence (plan, phase 7, step 35): tiers-novaseq's indexes (tbase05_all,tbase02_all,slim) queried again
# with the shipped coverage_zi and its interval, scored against true depth (abund_* in
# aai_summary.tsv). Needs the finished tiers-novaseq run; rebuilds no index.
[[ -d ../tiers-novaseq/results/mgnify_index ]] || { echo 'tiers-novaseq has not finished' >&2; exit 1; }
fmh abundance --aai_subset "$SUBSET" --prebuilt ../tiers-novaseq/results --prebuilt_names tbase05_all,tbase02_all,slim
