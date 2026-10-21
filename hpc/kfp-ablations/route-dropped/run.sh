#!/usr/bin/env bash
#SBATCH --job-name=kfp-route-dropped
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --export=ALL
# What routing (k = 8, min_kmers = 3) drops on hadza10M, searched with DIAMOND (--sensitive
# -F 15 --top 10) against every unit of the profile, with 1 and with 3 members per unit
# (plan, phase 13, routing), with workflows/route-dropped. Inputs: step 3b's unrouted
# sensitive_F15 arm in ../pileup-cost (its ref/, profile.tsv and pileup/assignments.parquet),
# hadza10M's reads, the standard tier's index and MERGE buckets. Reads are single-end: each mate routed on its own; 2 M
# dropped reads are aligned (READS=0: all). dropped.json: reads routed and dropped, by total
# k-mers.
# Read summary.tsv (a row per arm, m1 and m3, and per total-k-mer bin of the reads):
#   aligned, identity_* : dropped reads with an alignment, by best identity
#   best_on_*, high_on_*: the role of the unit hit (target, neighbour, profile only), at
#                         any identity and at >= 90%: high_on_target/neighbour are routing's
#                         misses; high_on_profile are reads of units the reference lacks
#   profile_units_high(_10), _median_q: the profile-only units those gather on
#   baseline_*          : how the unrouted run treated the same reads (aligned, on a unit, novel)
#   diamond_cpu_h, cpu_h_per_m_reads: the cost (makedb included), against step 3b's
#                         sensitive_F15 DIAMOND on all hadza10M pairs (1.79 CPU-h)
# and the log line "m3 vs m1": identity gained with 3 members.
source ../common.sh
# the workflow is $WF/route-dropped: a nextflow.config left here by its first layout would
# override its params (python = '')
[[ ! -f nextflow.config ]] || { echo "remove ./nextflow.config (stale; the workflow has its own)" >&2; exit 1; }
STD=../release-tiers/results-standard
base=../pileup-cost/results-dmnd-sensitive_F15/hadza10M
[[ -f $base/pileup/assignments.parquet ]] || { echo "no $base: run diamond-sensitivity first" >&2; exit 1; }
nextflow run "$WF/route-dropped" -c ../site.config -resume --outdir results \
    --ref "$(realpath "$base/ref")" --profile "$(realpath "$base/profile.tsv")" \
    --index "$(realpath "$STD/index")" --store "$(realpath "$STD")/members.bucket*.parquet" \
    --r1 "$(realpath ../pileup-cost/hadza10M_1.fastq.gz)" \
    --r2 "$(realpath ../pileup-cost/hadza10M_2.fastq.gz)" --reads "${READS:-2000000}"
"$REPO/.venv/bin/python" "$WF/route-dropped/dropped.py" summarise --dropped results/dropped \
    --arm m1=results/m1 --arm m3=results/m3 --assignments "$base/pileup/assignments.parquet" \
    --trace results/trace.tsv --out summary.tsv
