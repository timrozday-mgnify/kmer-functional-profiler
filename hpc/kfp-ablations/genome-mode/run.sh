#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-genome-mode
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Phase 11 (dev-phase11 branch): needs a checkout with `annotate-genomes` and `genomes`.
# Add --mgnify_full <dir> to also annotate with a built full MGnify index.
source ../common.sh
fmh genomes --aai_subset "$SUBSET"
