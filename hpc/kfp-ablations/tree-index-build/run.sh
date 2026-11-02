#!/usr/bin/env bash
#SBATCH --job-name=kfp-tree-index-build
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Tree index on the 1-in-100 subset (docs/tree-index-plan.md, Progress log: per-layer bars,
# gather, layer representatives and rough alignment), with workflows/tree-index/build.nf, then
# each tier benchmarked on BLOSUM62 strains.
#   base      S3, no layer (--n-1 2 --window-min 2 --min-kmers 32 --layers ''): t_base 0.01,
#             ~4 spread 11-mer minimizers, k 7 ladder (window 4) for leaves with < 32
#             informative 11-mers
#   slim      --mode slim: n_1 1, ~2 minimizers (--window-min 1) and the ladder for < 16;
#             ~100-130 GB at release (estimate)
#   standard  --mode standard: slim's base plus the rarest protein k-mer per window of 16 on
#             one composite per spanning-tree group (tree:0.57, k-mers of >= 25% of its
#             members), thinned to ~2.6 keys per group in proportion to sqrt of its members
#             (--layer-keys sqrt:2.6), k from the database's size (build.layer_k: k 8 here); a
#             representative's seed stands for its group; ~195-215 GB at release; ~0.73
#             recall at 80% aligned
#   large     --mode large: the layer at window 8, ~14 keys per group (sqrt:14); ~360-415 GB;
#             ~0.94 at 80%, ~0.86 at 70% aligned
#   prev-std  the previous standard (window 64, no thinning; ~235-260 GB, ~0.71 at 80%)
#   prev-large the previous large (base plus window 16; ~505-600 GB, ~0.82 at 70%)
# a leaf per distinct dense 8-mer set, the dense 8-mer topology, Fitch diffs on the spanning
# tree itself. PARTITION (~130 GB, one job) splits ~53 M proteins into parts of 1 M; BUILD runs
# per part and tier (~30 GB with a layer: 7.6-8.0 GB per 262 k proteins locally; 96 GB given); MERGE writes
# one index per tier. TARGETS/SIMULATE: 2,000 proteins, one strain each at 100-70% identity
# (BLOSUM62, depth lognormal mu 0.3 sigma 0.6); QUERY per tier with decoys (reads reversed)
# setting a 1% chance-share bar, then verify.py: each (read pair, leaf) with evidence aligned
# ungapped (BLOSUM62) against the leaf's representative, its own 1% decoy bar; SCORE, SUMMARY.
# Read:
#   results/bench.tsv                per tier: postings (x_base: over base), decoy bar and
#                                    seeds, candidates, recall at depth > 1 by identity
#                                    (verify_0.01_nb: aligned at the decoy bar, neighbour
#                                    rule; strategies.tsv also has ev_ and fam_ seeds alone)
#   results/runs/q-<tier>/verify.json  pairs aligned, seconds, pairs per second
#   results/index-<tier>/meta.json   postings by kind, bytes, positions, parts, build time
#   results/store-<tier>/            the leaf-ordered store (store.py, flat), read by verify.py
#                                    and tier2.py; ~1x the index's residues on disk
#   results/runs/s-<tier>/           recall.tsv (identity x depth), strategies.tsv
#   results/trace.tsv                time and peak RSS per task
# Re-submitting resumes; any change to the workflow's scripts or --part_proteins reruns from PARTITION.
source ../common.sh
[[ -f $SUBSET/members.parquet ]] || { echo "no $SUBSET/members.parquet" >&2; exit 1; }
tiers=(
    "base:--mode standard --layers '' --n-1 2 --window-min 2 --min-kmers 32"
    "slim:--mode slim"
    "standard:--mode standard"
    "large:--mode large"
    "prev-std:--mode standard --layers protein:{k}:64:rare --layer-keys ''"
    "prev-large:--mode standard --n-1 2 --window-min 2 --min-kmers 32 --layers protein:{k}:16:rare --layer-keys ''"
)
TIERS=$(IFS=';'; echo "${tiers[*]}")
# BUILD's memory grows with a part's proteins (~60 GB per 2 M with a layer, at the default
# 64 GB): parts of 1 M (~30 GB), 96 GB up front, +96 GB per retry. A component larger than a part
# (single linkage chains families: the largest outgrew 288 GB) is cut between MGnify90 clusters
# (partition.json: largest_component, split_components)
cat > build.config <<'CFG'
process {
    withName: '.*BUILD' { memory = { 96.GB * task.attempt } }
}
CFG
nextflow run "$WF/tree-index/build.nf" -c ../site.config -c build.config -resume --outdir results \
    --members "$(realpath "$SUBSET/members.parquet")" --pfam "$(realpath "$SUBSET/pfam.parquet")" \
    --tiers "$TIERS" --part_proteins 1000000
