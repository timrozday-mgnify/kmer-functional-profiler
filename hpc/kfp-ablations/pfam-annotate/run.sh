#!/usr/bin/env bash
#SBATCH --job-name=kfp-pfam-annotate
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# The benchmark genomes' Pfam annotations (proteins, hmmsearch --cut_ga, domains), once into
# the shared ../fmh-benchmark-data/pfam_annotations/; every fmh run with Pfam labels then reuses
# them instead of running hmmsearch again. Run before the next batch of fmh runs.
fmh '' --dbs_only --dbs pfam
