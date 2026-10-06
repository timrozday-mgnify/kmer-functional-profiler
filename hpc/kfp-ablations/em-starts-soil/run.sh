#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-em-starts-soil
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# EM convergence on deep, rich soil metagenomes, with no truth (plan, phase 7, step 30): each
# query cell fitted from 5 EM starts (the usual one and 4 random ones, --em-start) and the
# profiles compared in results-<run>/em_starts.tsv. Indexes: tiers-novaseq's slim, tbase02_all
# (standing in for the standard tier), tbase05_all and tbase10_all (whose largest components
# are fitted block-wise). Each run is profiled alone, not pooled, at 40 M, 100 M and all its
# pairs (2 x 151 bp, NovaSeq 6000):
#   SRR12659369  Harvard Forest LTER forest soil (PRJNA654925), 178,144,788 pairs, 24 GB
#   SRR11613032  Everglades Agricultural Area soil (PRJNA620814), 200,931,322 pairs, 26 GB
[[ -d ../tiers-novaseq/results/mgnify_index/tbase05_all/index ]] || { echo 'tiers-novaseq indexes missing' >&2; exit 1; }
# a config file, not --query_args: Nextflow reads a value starting with '--' as an option
printf 'params.query_args = "--aai"\n' > aai.config  # fits the zero-inflated EM (aai) too
for run in SRR12659369:178144788 SRR11613032:200931322; do
    nextflow run "$WF/mgnify-subset" -c ../site.config -c aai.config --query true --index false \
        --query_indexes '../tiers-novaseq/results/mgnify_index/{slim,tbase02_all,tbase05_all,tbase10_all}/index' \
        --query_run "${run%%:*}" --query_ladder "40000000,100000000,${run##*:}" \
        --query_starts 0,1,2,3,4 --query_draws 0 --outdir "results-${run%%:*}" -resume
done
