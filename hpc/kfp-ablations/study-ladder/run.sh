#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-study-ladder
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Phase 10 (dev branch): needs a checkout with workflows/study-index.
source ../common.sh
fmh study --aai_subset "$SUBSET"
