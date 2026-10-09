#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-unknown-novaseq
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# Unknown fraction (plan, phase 7, steps 14 and 36): slim and standard built on the 1-in-100
# subset and without the units 10/30/50% of the genomes hit, 10 NovaSeq metagenomes queried
# with --summary; read unknown_summary.tsv. Needs 00-dbs (Pfam); independent of the other runs.
fmh unknown --aai_subset "$SUBSET"
