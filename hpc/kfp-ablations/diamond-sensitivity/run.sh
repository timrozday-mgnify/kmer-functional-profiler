#!/usr/bin/env bash
#SBATCH --job-name=kfp-diamond-sensitivity
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --export=ALL
# DIAMOND sensitivity ablation (plan, phase 13, step 3b): --fast, default, --mid-sensitive and
# --sensitive (today's), each with -F 15 (today's) and -F 0. Each arm runs inside ../spike-protein
# (standard tier, truth) and ../pileup-cost (hadza10M, full standard tier, no truth), so their
# -resume caches give the samples, index, query and reference: only DIAMOND and what follows
# rerun (DIAMOND's options and -F are hashed inputs). Neither of those runs may be running.
# CliqueSNV, inStrain and LoFreq are off in every arm (freebayes and Floria kept), so the
# baseline arm reruns nothing but publishing. The head job holds compare.py's per-mate
# tables (hadza10M's assignments), hence 64 GB.
# Read compare.tsv (one row per benchmark, sample and arm, against sensitive_F15):
#   spike: read_precision_em, sibling_caught_em(_parent_target), l1_targets_em, f1 from
#     the arm's summary.tsv (truth);
#   every sample: mates aligned (arm / baseline), best-unit agreement per mate, em_abundance
#     Spearman and L1 per unit, freebayes record Jaccard, Floria multi-strain contigs;
#   cost: DIAMOND and pileup CPU-h, CALL's peak RSS (trace.tsv).
# Then decision.tsv: the rule fixed before the run (plan, step 3b), among the -F 15 arms only
# (frameshifts are a wanted output; the -F 0 arms are for cost): read precision >= baseline
# - 0.005, siblings caught (EM) >= baseline - 0.02, em_abundance Spearman >= 0.99 and
# freebayes Jaccard >= 0.98 in every sample; the cheapest passing arm by DIAMOND CPU-h.
source ../common.sh
# their run.sh writes these: they must have run (their caches are this run's input)
[[ -f ../spike-protein/spike.config && -f ../pileup-cost/cost.config ]] ||
    { echo 'run spike-protein and pileup-cost first' >&2; exit 1; }
here=$PWD
ARMS=${ARMS:-sensitive_F15 mid-sensitive_F15 default_F15 fast_F15 sensitive_F0 mid-sensitive_F0 default_F0 fast_F0}
for arm in $ARMS; do
    mode=${arm%_F*} f=${arm##*_F}
    # a config, not the command line: Nextflow reads a value starting with '--' as an option
    printf 'params {\n    diamond_args = "%s"\n    frameshift = %s\n    cliquesnv = false\n    instrain = false\n    lofreq = false\n}\n' \
        "$([[ $mode == default ]] || echo "--$mode")" "$f" > "$here/arm-$arm.config"
done
for arm in $ARMS; do
    (cd ../spike-protein && nextflow run "$WF/spike-benchmark" -c ../site.config -c spike.config \
        -c "$here/arm-$arm.config" -resume \
        --seed_fasta seeds.faa --families 2000 --seeds 10 \
        --siblings_args '--siblings 0.1' --background "$PWD/bg_1.fastq.gz,$PWD/bg_2.fastq.gz" \
        --index_args "--k 11 --t-base 0.02 --n-min 16 --oversample 1" \
        --aai_model "$WF/fmh-benchmark/ablations/aai_models/step34q.json" \
        --query_stats true --outdir "results-dmnd-$arm")
done
STD=../release-tiers/results-standard
printf 'sample,r1,r2\nhadza10M,hadza10M_1.fastq.gz,hadza10M_2.fastq.gz\n' > ../pileup-cost/samples-dmnd.csv
for arm in $ARMS; do
    (cd ../pileup-cost && nextflow run "$WF/pileup" -c ../site.config -c cost.config \
        -c "$here/arm-$arm.config" -resume --samples samples-dmnd.csv \
        --index "$STD/index" --store "$STD/members.bucket*.parquet" --outdir "results-dmnd-$arm")
done
"$REPO/.venv/bin/python" compare.py --base sensitive_F15 --arms $ARMS \
    --bench spike=../spike-protein --bench real=../pileup-cost
