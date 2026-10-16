#!/usr/bin/env bash
#SBATCH --job-name=kfp-ablation-g16-novaseq
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# max_groups 16 on the slim and standard tiers (plan, phase 7, step 43): each built at 64 (control)
# and 16 on the 1-in-100 subset, 10 NovaSeq metagenomes; read summary.tsv and aai_summary.tsv
# against the 64 controls. tbase05_all_g16 again, against promiscuity-novaseq's, for the
# zero-inflated fit's 10,000-step cap. Needs 00-dbs (Pfam); before release-tiers (builds at 16).
fmh g16 --aai_subset "$SUBSET"
