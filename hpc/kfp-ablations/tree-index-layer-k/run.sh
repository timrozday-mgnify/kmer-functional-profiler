#!/usr/bin/env bash
#SBATCH --job-name=kfp-tree-index-layer-k
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Tree index, the protein layer's k at the 1-in-100 subset's size (docs/tree-index-plan.md,
# Progress log: tree-index-build results). There, layer k 8 (build.layer_k's log-size rule) lost
# recall below ~85% against the quarter sample's k 6: standard 0.47 / 0.23 at 80 / 70% aligned
# (quarter 0.73 / 0.58). Same subset, partition and benchmark as tree-index-build (2,000 BLOSUM62
# strains at 100-70%), tiers differing in the layer only:
#   standard    k 8, ~2.6 keys per group (sqrt:2.6): the control
#   std-k7      k 7
#   std-k6      k 6, the quarter's k (chance seeds ~400x k 8's: verify.py's pairs and QUERY's
#               memory are the risk)
#   std-keys2x  k 8, twice the keys (sqrt:5.2): more seeds per group at k 8's chance rate
#   large       k 8, window 8, ~14 keys per group: the control
#   large-k7    k 7
# A tier whose QUERY fails twice is skipped (errorStrategy ignore), so the rest still score.
# Read:
#   results/bench.tsv                 per tier: postings, bytes, decoy bar, recall by identity,
#                                     tier 2's read recall and precision, Pfam precision (over
#                                     calls whose source has a Pfam hit) and pfam_known_share
#   results/runs/q-<tier>/verify.json pairs aligned, seconds; tier2/tier2.json (routing, pairs)
#   results/runs/s-<tier>/            reads.tsv and detection.tsv, tier 1 and tier2/
#   results/trace.tsv                 time and peak RSS per task (QUERY's for k 6 especially)
# Re-submitting resumes.
source ../common.sh
[[ -f $SUBSET/members.parquet ]] || { echo "no $SUBSET/members.parquet" >&2; exit 1; }
tiers=(
    "standard:--mode standard"
    "std-k7:--mode standard --layers protein:7:16:rare"
    "std-k6:--mode standard --layers protein:6:16:rare"
    "std-keys2x:--mode standard --layer-keys sqrt:5.2"
    "large:--mode large"
    "large-k7:--mode large --layers protein:7:8:rare"
)
TIERS=$(IFS=';'; echo "${tiers[*]}")
# BUILD as tree-index-build's (parts of 1 M; 85-91 GB on the giant component's parts there)
cat > build.config <<'CFG'
process {
    withName: '.*BUILD' { memory = { 96.GB * task.attempt } }
    withName: '.*QUERY' { errorStrategy = { task.attempt < 2 ? 'retry' : 'ignore' } }
}
CFG
nextflow run "$WF/tree-index/build.nf" -c ../site.config -c build.config -resume --outdir results \
    --members "$(realpath "$SUBSET/members.parquet")" --pfam "$(realpath "$SUBSET/pfam.parquet")" \
    --tiers "$TIERS" --part_proteins 1000000
