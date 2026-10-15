#!/usr/bin/env bash
#SBATCH --job-name=kfp-spike-nucleotide
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Spike-in benchmark, nucleotide indexes (plan, phase 14, step 5 and its gate): CARD's
# protein-homolog nucleotide models as targets, strains at 100/98/95/90% identity, 10% with
# a withheld sibling; 10 samples; bwa-mem2 follow-up. Arms (results-<arm>/summary.tsv, rows
# kfp, kfp_aa, kma):
#   k21         k 21, t_base 0.02, n_min 16, with KMA and kfp on the translated targets: the
#               gate (F1 and L1 within 0.05 of KMA's, ani |median bias| <= 0.01 at 0.9-1.0)
#   k21_bg      as k21, spiked into the first 10 M pairs of ERR7738575: real false positives
#               (gut ARGs count against both tools) and cost against KMA's (trace.tsv)
#   k17 k25 k31 k21_t10 k21_dense   tuning (k, t_base; plan: "re-tuned on the nucleotide
#               benchmark rather than carried over"), kfp only
# CARD is free for academic use (https://card.mcmaster.ca/about); the download is yours.
source ../common.sh
if [[ ! -s card.fna ]]; then
    curl -fsSL https://card.mcmaster.ca/latest/data -o card.tar.bz2
    mkdir -p card && tar -xjf card.tar.bz2 -C card
    # ids without spaces (the first word is kept by every tool)
    sed 's/ /_/g' "$(find card -name nucleotide_fasta_protein_homolog_model.fasta | head -n 1)" > card.fna
fi
cat > spike.config <<CFG
process {
    withName: 'SIMULATE' { memory = 16.GB; time = 8.h }
    withName: 'PILEUP:BWAMEM2|KMA' { cpus = 16; memory = { 64.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: 'PILEUP:CALL' { memory = { 64.GB * task.attempt } }
}
CFG
common=(--molecule nucleotide --targets card.fna --seeds 10 --siblings_args '--siblings 0.1'
        --query_stats true -c ../site.config -c spike.config -resume)
for arm in "k21:--k 21 --t-base 0.02 --n-min 16 --oversample 1:true:" \
           "k21_bg:--k 21 --t-base 0.02 --n-min 16 --oversample 1:true:bg" \
           "k17:--k 17 --t-base 0.02 --n-min 16 --oversample 1:false:" \
           "k25:--k 25 --t-base 0.02 --n-min 16 --oversample 1:false:" \
           "k31:--k 31 --t-base 0.02 --n-min 16 --oversample 1:false:" \
           "k21_t10:--k 21 --t-base 0.1 --n-min 16 --oversample 1:false:" \
           "k21_dense:--k 21 --t-base 1 --n-min 0:false:"; do
    IFS=: read -r name args others bg <<< "$arm"
    extra=()
    if [[ -n $bg ]]; then
        fetch_pairs bg 10000000 ERR7738575
        extra=(--background "$PWD/bg_1.fastq.gz,$PWD/bg_2.fastq.gz")
    fi
    nextflow run "$WF/spike-benchmark" "${common[@]}" ${extra[@]+"${extra[@]}"} \
        --index_args "$args" --kma "$others" --aa_arm "$others" \
        --aa_index_args '--k 11 --t-base 0.02 --n-min 16 --oversample 1' --outdir "results-$name"
done
