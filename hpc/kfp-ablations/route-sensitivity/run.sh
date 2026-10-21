#!/usr/bin/env bash
#SBATCH --job-name=kfp-route-sensitivity
#SBATCH --output=%x.%j.log
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --export=ALL
# Routing's k and min_kmers against unrouted DIAMOND (plan, phase 13, routing). No DIAMOND
# run: the baseline is step 3b's finished sensitive_F15 arm (unrouted, --sensitive -F 15) in
# ../spike-protein (10 samples, truth) and ../pileup-cost (hadza10M); sweep.py routes every
# sample's pairs at component granularity for each k and thresholds at each min_kmers.
# Read decision.tsv (per k, min_kmers: the worst sample's unit_recall, best_recall,
# strain_recall and sibling_recall; mean routed_share, searched_share and novel_routed;
# passes: the rule in sweep.py, fixed before the run) and sweep.tsv per sample; hist.*.tsv:
# pairs by the most k-mers shared with one component. Chance hits grow as the reference's
# distinct k-mers over 20^k: at k = 7 against hadza10M's ~20 M residues ~1.5% of a read's
# k-mers hit by chance, so small k needs a larger min_kmers. Memory: hadza10M's assignments (23.6 M rows) and
# per-pair hits at k = 7.
source ../common.sh
PY="$REPO/.venv/bin/python"
K=${K:-7 8 9 10 11}
MIN=${MIN:-1 2 3 5 8}
mkdir -p out

# a spike task's work directory, from the arm's trace (its inputs are not published)
spike=../spike-protein
base=$spike/results-dmnd-sensitive_F15
[[ -f $base/trace.tsv ]] || { echo "no $base: run diamond-sensitivity first" >&2; exit 1; }
column() {  # trace.tsv rows of finished tasks: name, hash
    awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
        $c["status"] == "COMPLETED" || $c["status"] == "CACHED" { print $c["name"] "\t" $c["hash"] }' \
        "$base/trace.tsv"
}
workdir() {  # a task's exact name -> its work directory
    local hash
    hash=$(column | awk -F '\t' -v n="$1" '$1 == n { h = $2 } END { print h }')
    [[ -n $hash ]] || { echo "no finished task '$1' in $base/trace.tsv" >&2; exit 1; }
    ls -d "$spike"/work/"$hash"*
}
index_task=$(column | awk -F '\t' 'index($1, "INDEX_TARGETS") == 1 { print $1; exit }')
index=$(find "$(workdir "$index_task")"/ -maxdepth 1 -type d -name 'index_*' | head -n 1)
for s in $(cd "$base" && ls -d seed*/pileup | cut -d/ -f1 | sort -V); do
    w=$(workdir "SIMULATE (seed ${s#seed})")
    [[ -s out/sweep.spike.$s.tsv ]] ||
        "$PY" sweep.py --bench spike --sample "$s" --ref "$base/$s/ref" \
            --assignments "$base/$s/pileup/assignments.parquet" --index "$index" \
            --r1 "$w/reads_1.fq.gz" --r2 "$w/reads_2.fq.gz" --truth "$w/truth.tsv" \
            --k $K --min-kmers $MIN --out out
done
real=../pileup-cost/results-dmnd-sensitive_F15/hadza10M
[[ -s out/sweep.real.hadza10M.tsv ]] ||
    "$PY" sweep.py --bench real --sample hadza10M --ref "$real/ref" \
        --assignments "$real/pileup/assignments.parquet" \
        --index ../release-tiers/results-standard/index \
        --r1 ../pileup-cost/hadza10M_1.fastq.gz --r2 ../pileup-cost/hadza10M_2.fastq.gz \
        --k $K --min-kmers $MIN --out out
"$PY" sweep.py decide out
