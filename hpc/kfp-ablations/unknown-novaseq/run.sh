#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-unknown-novaseq
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# Unknown fraction (plan, phase 7, steps 14, 36 and 40): slim and standard built on the 1-in-100
# subset and without the units 10/30/50% of the sampled genomes hit, 10 NovaSeq metagenomes
# queried with --summary; read census_containment in unknown_summary.tsv. Needs 00-dbs (Pfam);
# independent of the other runs. Step 40 resubmits it: -resume reuses the DIAMOND passes.
fmh unknown --aai_subset "$SUBSET"
