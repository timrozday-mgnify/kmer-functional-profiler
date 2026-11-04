#!/usr/bin/env bash
#SBATCH --job-name=kfp-tree-index-benchmark
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Tree index, phase 4 (docs/tree-index-plan.md): tiers already built, against today's tier 1 on
# the same proteins and the same strains, scored without fixed units, with
# workflows/tree-index/build.nf --built. 2,000 BLOSUM62 strains at 100-70%, no background.
#   TARGETS..SCORE  as tree-index-build, neighbour rule at delta 0.05, 0.1 and 0.15
#   FLAT      today's kfp index per tier (slim, standard, large), kfp query --aai
#   REPORTED  tree tier 1, tree tier 2 and flat: units expanded to members
#   ALIGN     DIAMOND, strains against reported: identities independent of the methods
#   NEAREST   curves, earth mover's distance, roll-ups (MGnify90, Pfam), AAI
#   GATE      each tree tier and its tier 2 against the flat tier of its name
# BUILT=<a build run's results> (default: tree-index-build's) and TIERS='standard:;large:' pick
# the index. KOfam: pass its table with nearest.py score --rollup kofam=<parquet>.
# Read:
#   results/gate.tsv                    per run and measure: tier 1, tree at tier 1's precision, pass
#   results/nearest/<run>/score.json    read-outs at 99, 95, 90, 80%; EMD; roll-ups; cost
#   results/nearest/<run>/curves.tsv    nearest-match curves, overall, by identity and depth
#   results/nearest/<run>/aai.tsv       estimated AAI against true identity
#   results/nearest/tree-<tier>/neighbours.tsv  spikes with no own hit reached by the rule
#   results/trace.tsv                   time and peak RSS per task
# Re-submitting resumes.
source ../common.sh
BUILT=${BUILT:-../tree-index-build/results}
TIERS=${TIERS:-'slim:;standard:;large:'}
[[ -f $SUBSET/members.parquet ]] || { echo "no $SUBSET/members.parquet" >&2; exit 1; }
for t in ${TIERS//;/ }; do
    [[ -d $BUILT/index-${t%%:*} && -d $BUILT/store-${t%%:*} ]] || { echo "no $BUILT/{index,store}-${t%%:*}" >&2; exit 1; }
done
nextflow run "$WF/tree-index/build.nf" -c ../site.config -resume --outdir results \
    --built "$(realpath "$BUILT")" --tiers "$TIERS" --score_deltas '0.05 0.1 0.15' \
    --members "$(realpath "$SUBSET/members.parquet")" --pfam "$(realpath "$SUBSET/pfam.parquet")"
