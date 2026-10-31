#!/usr/bin/env bash
#SBATCH --job-name=kfp-tree-index-build
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Tree index on the 1-in-100 subset (docs/tree-index-plan.md, Progress log: tiers by identity),
# with workflows/tree-index/build.nf, then each tier benchmarked on BLOSUM62 strains.
#   base      S3, no layer (--mode standard --layers ''): t_base 0.01, ~4 spread 11-mer
#             minimizers, k 7 ladder (window 4) for leaves with < 32 informative 11-mers
#   slim      --mode slim: as base with n_1 1, ~2 minimizers (--window-min 1) and the ladder
#             for < 16; ~100-130 GB at release (estimate)
#   standard  --mode standard: base plus the rarest protein k-mer per window of 8 on every leaf,
#             k from the database's size (build.layer_k: k 8 here, ~32 chance bits); recall to
#             ~80% identity
#   large     --mode large: the same layer at window 4; recall to ~70%
#   std-w16   the layer at window 16, between base and standard
# a leaf per distinct dense 8-mer set, the dense 8-mer topology, Fitch diffs on the spanning
# tree itself. PARTITION (~130 GB, one job) splits ~53 M proteins into parts of 2 M; BUILD runs
# per part and tier (~60 GB with a layer: 7.6-8.0 GB per 262 k proteins locally); MERGE writes
# one index per tier. TARGETS/SIMULATE: 2,000 proteins, one strain each at 100-70% identity
# (BLOSUM62, depth lognormal mu 0.3 sigma 0.6); QUERY per tier with decoys (reads reversed)
# setting a 1% chance-share bar; SCORE and SUMMARY.
# Read:
#   results/bench.tsv                per tier: postings (x_base: over base), decoy bar and
#                                    seeds, candidates, recall at depth > 1 by identity
#                                    (ev_0.01_nb: seeds at the decoy bar, neighbour rule)
#   results/index-<tier>/meta.json   postings by kind, bytes, positions, parts, build time
#   results/runs/s-<tier>/           recall.tsv (identity x depth), strategies.tsv
#   results/trace.tsv                time and peak RSS per task
# Re-submitting resumes: PARTITION is reused; BUILD and MERGE rerun when build.py changes.
source ../common.sh
[[ -f $SUBSET/members.parquet ]] || { echo "no $SUBSET/members.parquet" >&2; exit 1; }
tiers=(
    "base:--mode standard --layers ''"
    "slim:--mode slim"
    "standard:--mode standard"
    "large:--mode large"
    "std-w16:--mode standard --layers protein:{k}:16:rare"
)
TIERS=$(IFS=';'; echo "${tiers[*]}")
nextflow run "$WF/tree-index/build.nf" -c ../site.config -resume --outdir results \
    --members "$(realpath "$SUBSET/members.parquet")" --tiers "$TIERS"
