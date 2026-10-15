#!/usr/bin/env bash
#SBATCH --job-name=kfp-pileup-cost
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# The pileup follow-up on real reads against the full standard tier (plan, phase 13: steps
# 0, 3 and 6; open questions: DIAMOND's wall time against pass 1 at 10 M and 100 M pairs).
# Samples: the first 10 M pairs of ERR7738575 (Hadza gut) and the first 100 M of ERR7738575 +
# ERR7746321 pooled (phase 6's pair). Needs the finished release-tiers run (its
# results-standard/index, step34q attached, and the MERGE buckets as the sequence store).
# Read:
#   trace.tsv           PILEUP:DIAMOND against PILEUP:QUERY realtime: step 3's rule (> ~5x
#                       at 10 M pairs: add the in-memory dense pre-filter)
#   <sample>/stats.json hits_per_pair: step 0's histogram (replaces the plan's estimates)
#   <sample>/aai_check.tsv  step 6 on real data: aai against the reads' alignment identity
#                       (flagged rows are phase-7 aai findings)
#   <sample>/variants/  step 6c on real data: freebayes calls, csq consequences, Floria and
#                       CliqueSNV haplotypes, inStrain microdiversity against the consensus (6b)
#   results-noF/trace.tsv  PILEUP:DIAMOND without -F (step 6b's cost: against trace.tsv)
# Buckets written before phase 13, step 1 (10,000-row groups) make fetch read far more.
source ../common.sh
STD=../release-tiers/results-standard
[[ -f $STD/index/meta.json ]] || { echo 'release-tiers has not finished (no results-standard/index)' >&2; exit 1; }
fetch_pairs hadza10M 10000000 ERR7738575
fetch_pairs pooled100M 100000000 ERR7738575 ERR7746321
printf 'sample,r1,r2\nhadza10M,hadza10M_1.fastq.gz,hadza10M_2.fastq.gz\npooled100M,pooled100M_1.fastq.gz,pooled100M_2.fastq.gz\n' > samples.csv
cat > cost.config <<CFG
process {
    withName: 'PILEUP:QUERY' { cpus = 8; memory = { 64.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: 'PILEUP:REFERENCE' { memory = { 64.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: 'PILEUP:DIAMOND' { cpus = 32; memory = { 128.GB * task.attempt }; time = { 72.h * task.attempt } }
    withName: 'PILEUP:CALL' { memory = { 128.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: 'PILEUP:MAP_CONSENSUS' { cpus = 16; memory = { 32.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: 'PILEUP:(HANDOFF|BAM|FREEBAYES|CONSEQUENCES)' { memory = { 16.GB * task.attempt }; time = { 48.h * task.attempt } }
    withName: 'PILEUP:(FLORIA|INSTRAIN|CLIQUESNV)' { cpus = 16; memory = { 64.GB * task.attempt }; time = { 48.h * task.attempt } }
}
CFG
nextflow run "$WF/pileup" -c ../site.config -c cost.config -resume --samples samples.csv \
    --index "$STD/index" --store "$STD/members.bucket*.parquet" --query_stats true --outdir results
# step 6b's cost of frameshift-aware alignment: the 10 M sample again with -F off
printf 'sample,r1,r2\nhadza10M,hadza10M_1.fastq.gz,hadza10M_2.fastq.gz\n' > samples-noF.csv
nextflow run "$WF/pileup" -c ../site.config -c cost.config -resume --samples samples-noF.csv \
    --index "$STD/index" --store "$STD/members.bucket*.parquet" --frameshift 0 \
    --lofreq false --floria false --instrain false --cliquesnv false \
    --outdir results-noF
for s in hadza10M pooled100M; do
    "$REPO/.venv/bin/kfp-pileup" aai-check "results/$s/pileup/units.tsv" --out "results/$s/aai_check.tsv"
done
