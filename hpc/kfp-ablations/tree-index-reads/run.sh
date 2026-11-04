#!/usr/bin/env bash
#SBATCH --job-name=kfp-tree-index-reads
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Tree index, tier 1 and tier 2 on a real sample (docs/tree-index-plan.md, phase 4b): hadza10M
# (pileup-cost's 10 M read pairs of ERR7738575) through tiers already built, with
# workflows/tree-index/build.nf --reads. No build and no truth: cost at real-sample scale.
#   TIER1   per tier: query.py --decoy (reads reversed for the 1% chance-share bar), verify.py
#   TIER2   per tier and memory budget (16 and 4 GB): tier2.py (dense routing, two-stage
#           alignment, EM over ties, the hand-off in kfp-pileup reference's format)
#   CHUNKS  per tier: the two budgets' calls compared; at most 1% of called reads may differ
# BUILT=<a build run's results> (default: tree-index-build's, the 1-in-100 subset; later
# tree-index-release's) and TIERS='standard:;large:' pick the index. Compare against the flat
# table on the same reads: ../pileup-fdr/results/hadza10M/ and its trace.tsv.
# Read:
#   results/q-<tier>/query.json, verify.json   tier 1: seeds, decoy bar, pairs aligned, seconds
#   results/t2-<tier>-<gb>gb/tier2.json        tier 2: reads routed, pairs, chunks, batches,
#                                              seconds per step; reads.parquet (calls);
#                                              handoff/ (ref.fa, refs.tsv, assignments.parquet,
#                                              reads_{1,2}.fq.gz) for kfp-pileup's DIAMOND
#   results/chunks-<tier>.json                 identical across budgets (exact, to rounding), share of reads differing
#   results/trace.tsv                          time and peak RSS per task
# Re-submitting resumes.
source ../common.sh
BUILT=${BUILT:-../tree-index-build/results}
TIERS=${TIERS:-'slim:;standard:;large:'}
for t in ${TIERS//;/ }; do
    [[ -d $BUILT/index-${t%%:*} && -d $BUILT/store-${t%%:*} ]] || { echo "no $BUILT/{index,store}-${t%%:*}" >&2; exit 1; }
done
fetch_pairs hadza10M 10000000 ERR7738575
nextflow run "$WF/tree-index/build.nf" -c ../site.config -resume --outdir results \
    --reads "$(realpath hadza10M_1.fastq.gz),$(realpath hadza10M_2.fastq.gz)" \
    --built "$(realpath "$BUILT")" --tiers "$TIERS" --tier2_memory 16,4
