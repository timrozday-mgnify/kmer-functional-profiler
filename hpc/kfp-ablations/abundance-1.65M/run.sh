#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-abundance-1.65M
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# Abundance under divergence (plan, phase 7, step 35): floor-1.65M's indexes (tbase05,tbase02_all,dense02) queried again
# with the shipped coverage_zi and its interval, scored against true depth (abund_* in
# aai_summary.tsv). Needs the finished floor-1.65M run; rebuilds no index.
[[ -d ../floor-1.65M/results/mgnify_index ]] || { echo 'floor-1.65M has not finished' >&2; exit 1; }
fmh abundance --aai_subset "$SUBSET" --prebuilt ../floor-1.65M/results --prebuilt_names tbase05,tbase02_all,dense02 --n_reads 1650000
