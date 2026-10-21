#!/usr/bin/env bash
#SBATCH --job-name=kfp-pileup-fdr
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# The pileup on hadza10M with presence-fdr's targets (plan, phase 13: presence at a false
# discovery rate; route-dropped's "Next"): query at the defaults (--background-share 0.1),
# targets at present_q <= 0.05, routing at (8, 3), DIAMOND --sensitive -F 15. Its own
# directory, so nothing is reused from pileup-cost's cache (whose QUERY predates presence-fdr);
# pileup-cost's reads are read in place. MAX_TARGETS=N caps targets (the most hit first) if
# REFERENCE, DIAMOND or CALL is out of reach uncapped.
# Read:
#   results/hadza10M/profile.tsv  units at present_q <= 0.05 (simulated: ~2.2 M) and their
#                                 share of abundance_zi (before: 9,931, 1.1%)
#   results/hadza10M/route.json   reads routed and dropped (before: 22% of reads routed)
#   trace.tsv                     REFERENCE, ROUTE, DIAMOND and CALL CPU-h and peak RSS against
#                                 pileup-cost's hadza10M (DIAMOND 1.79 CPU-h unrouted, 32 k units)
#   results/hadza10M/call_stats.json, pileup/units.tsv  CALL's stages; units with reads
# Then route-dropped on it: cd ../route-dropped && BASE=../pileup-fdr/results/hadza10M OUT=results-fdr sbatch run.sh
source ../common.sh
STD=../release-tiers/results-standard
[[ -f $STD/index/meta.json ]] || { echo 'release-tiers has not finished (no results-standard/index)' >&2; exit 1; }
[[ -s ../pileup-cost/hadza10M_1.fastq.gz ]] || { echo 'no ../pileup-cost/hadza10M_*.fastq.gz: run pileup-cost first' >&2; exit 1; }
printf 'sample,r1,r2\nhadza10M,%s,%s\n' "$(realpath ../pileup-cost/hadza10M_1.fastq.gz)" \
    "$(realpath ../pileup-cost/hadza10M_2.fastq.gz)" > samples.csv
cat > cost.config <<CFG
process {
    resourceLimits = [memory: 512.GB]
    withName: 'PILEUP:QUERY' { cpus = 8; memory = { 64.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: 'PILEUP:REFERENCE' { memory = { 128.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: 'PILEUP:ROUTE' { cpus = 8; memory = { 64.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: 'PILEUP:DIAMOND' { cpus = 32; memory = { 128.GB * task.attempt }; time = { 72.h * task.attempt } }
    withName: 'PILEUP:CALL' { memory = { 256.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: 'PILEUP:MAP_CONSENSUS' { cpus = 16; memory = { 32.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: 'PILEUP:(HANDOFF|BAM|FREEBAYES|CONSEQUENCES)' { memory = { 16.GB * task.attempt }; time = { 48.h * task.attempt } }
    withName: 'PILEUP:(FLORIA|INSTRAIN|CLIQUESNV)' { cpus = 16; memory = { 64.GB * task.attempt }; time = { 48.h * task.attempt } }
}
CFG
# CliqueSNV, inStrain and LoFreq off, as step 3b's arms: the cost read here is the reference's
nextflow run "$WF/pileup" -c ../site.config -c cost.config -resume --samples samples.csv \
    --index "$STD/index" --store "$STD/members.bucket*.parquet" --query_stats true \
    --max_targets "${MAX_TARGETS:-0}" --cliquesnv false --instrain false --lofreq false \
    --outdir results
