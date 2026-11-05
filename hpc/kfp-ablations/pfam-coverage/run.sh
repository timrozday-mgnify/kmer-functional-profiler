#!/usr/bin/env bash
#SBATCH --job-name=kfp-pfam-coverage
#SBATCH --output=%x.%j.log
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# How many MGnify90 clusters have a Pfam hit on any member (plan, phase 6): one PFAM_COVERAGE
# job over mgy_cluster_seqs, mgy_proteins_pfam and mgy_clusters (~160 GB read, no extraction).
#   SAMPLE=100 sbatch run.sh     an estimate from 1 in 100 clusters (default 1: all of them)
#   RELEASE=<mirror> sbatch run.sh     reads a local mirror instead of the FTP site
# Read:
#   results/pfam_coverage.tsv    all, singleton and non-singleton clusters: clusters, with_pfam,
#                                fraction. At SAMPLE=1 clusters should be 1,661,142,385 in all.
#   results/trace.tsv            time and peak RSS
# Re-submitting resumes.
source ../common.sh
nextflow run "$WF/mgnify-subset" -c ../site.config ${RELEASE:+--release "$RELEASE"} \
    --pfam_coverage true --index false --sample "${SAMPLE:-1}" --outdir results -resume
