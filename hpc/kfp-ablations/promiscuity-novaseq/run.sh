#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-promiscuity-novaseq
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# Promiscuity cut (plan, phase 7, step 40): tbase05_all at max_groups 64/128/16/8/4 and tbase10_all
# at 64/8 built on the 1-in-100 subset, 10 NovaSeq metagenomes; read summary.tsv and
# aai_summary.tsv against the 64 controls. Needs 00-dbs (Pfam); independent of the other runs.
fmh promiscuity --aai_subset "$SUBSET"
