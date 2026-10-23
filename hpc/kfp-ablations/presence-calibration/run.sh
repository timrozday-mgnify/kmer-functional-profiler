#!/usr/bin/env bash
#SBATCH --job-name=kfp-presence-calibration
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# query --background-share calibrated at full-index scale (plan, phase 13 / phase 7: presence
# at a false discovery rate; the calibration design of 2026-10-21), with
# workflows/presence-calibration. 10 seeds x 2,000 units of the release standard tier, strains
# at 100/98/95/90% and ~0.05-20x spiked into hadza10M (pileup-cost's reads); hadza10M itself
# (unspiked) and reversed (null); shares 0.01, 0.03, 0.1, 0.3 and inf: 60 queries.
# Read:
#   results/decision.tsv  per share at q = 0.05: passes (every seed false_share <= 0.05 and
#                         null_called <= 0.01 x unspiked_called), recall; the chosen share
#   results/summary.tsv   per share and q (0.01, 0.05, 0.1): worst seed's false_share, null
#                         calls, recall by depth, identity and uniqueness bins
#   results/scores.tsv    per seed; results/<sample>/share<s>/stats.json, trace.tsv: query cost
# The rule (above) is params q, max_false and max_null: fix it before reading the results.
# Hit-level re-run (2026-10-23: recall before presence, for a routing threshold):
#   OUT=results-aai AAI=true sbatch run.sh
# query --aai lists every hit unit (~8 M rows, ~2-3 GB a profile; 60 profiles, copied out of
# work/ into $OUT). Read summary.tsv's recall_hit1, recall_hit2 and recall_gathered beside
# recall, recall_hit1_depth_* and _identity_*, and unspiked_hit*, null_hit* (reference size,
# chance hits); $OUT/seed<N>/truth.tsv is now published.
# Since 2026-10-23 PICK keeps members >= 0.9 of the representative's length (--min-member-
# frac), so runs before then (results, results-aai) are not comparable with later ones.
source ../common.sh
STD=../release-tiers/results-standard
[[ -f $STD/index/meta.json ]] || { echo 'release-tiers has not finished (no results-standard/index)' >&2; exit 1; }
[[ -s ../pileup-cost/hadza10M_1.fastq.gz ]] || { echo 'no ../pileup-cost/hadza10M_*.fastq.gz: run pileup-cost first' >&2; exit 1; }
nextflow run "$WF/presence-calibration" -c ../site.config -resume --outdir "${OUT:-results}" \
    --aai "${AAI:-false}" \
    --index "$(realpath "$STD/index")" --store "$(realpath "$STD")/members.bucket*.parquet" \
    --background "$(realpath ../pileup-cost/hadza10M_1.fastq.gz),$(realpath ../pileup-cost/hadza10M_2.fastq.gz)"
