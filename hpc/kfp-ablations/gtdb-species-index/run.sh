#!/usr/bin/env bash
#SBATCH --job-name=kfp-gtdb-species-index
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Phase 11, step 11 (dev-species branch): a GTDB r226 species index, genomes from NCBI,
# genes by pyrodigal. INDEX=<the built index dir profiles use> is required. Defaults to 1000
# sampled species (~0.5 TB less download than the full release); --max_species '' for all.
source ../common.sh
: "${INDEX:?INDEX=<built index dir>}"
nextflow run "$WF/gtdb-species-index" -c ../site.config --index "$INDEX" \
    --max_species 1000 --outdir results -resume "$@"
