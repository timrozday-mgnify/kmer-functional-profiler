#!/usr/bin/env bash
#SBATCH --job-name=kfp-pfam-reannotate
#SBATCH --output=%x.%j.log
#SBATCH --time=14-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# MGnify proteins against Pfam-A 38.2 with hmmsearch, keeping hits below GA (plan, phase 6;
# workflows/pfam-reannotate/README.md has the expected times).
#   sbatch run.sh                    1 in 100 cluster representatives (16.6 M), each searched
#                                    forwards and reversed, the null arm, keeping posterior
#                                    probabilities (~300 core-hours)
#   SAMPLE=1 sbatch run.sh           all representatives, forwards only, alignments without
#                                    posterior probabilities (~14,900 core-hours)
#   RELEASE=<mirror> sbatch run.sh   reads a local mirror instead of the FTP site
#   HMM=<Pfam-A.hmm[.gz]> sbatch run.sh    a local Pfam-A file instead of downloading it
# Use a new directory or clear results/ when changing SAMPLE.
# Read:
#   results/summary.tsv     per arm: proteins, residues, domains, domains_above_ga,
#                           proteins_with_hit, proteins_above_ga
#   results/compare.tsv     with the null arm: domains per 5-bit score bin, real and reversed
#   results/hits/*.parquet  one row per domain, with above_ga and the alignment (cigar, pp)
#   results/ga.tsv          each model's GA thresholds
#   results/trace.tsv       time and peak RSS per task
# Re-submitting resumes. work/ holds FASTA and HMMER's output until deleted (see the README).
source ../common.sh
sample=${SAMPLE:-100}
sampled=$([[ $sample == 1 ]] && echo false || echo true)  # the null arm and pp: 1-in-N runs only
nextflow run "$WF/pfam-reannotate" -c ../site.config ${RELEASE:+--release "$RELEASE"} \
    ${HMM:+--hmm "$HMM"} --targets reps --sample "$sample" \
    --reversed "$sampled" --pp "$sampled" \
    --outdir results -resume
