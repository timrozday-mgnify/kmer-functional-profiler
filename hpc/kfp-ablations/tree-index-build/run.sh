#!/usr/bin/env bash
#SBATCH --job-name=kfp-tree-index-build
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Tree index on the 1-in-100 subset (docs/tree-index-plan.md, phase 6b: reduced amino-acid
# alphabets), with workflows/tree-index/build.nf, then each tier benchmarked on BLOSUM62 strains.
#   standard  build.py --mode standard: t_base 0.01, no floors, ~4 spread 11-mer minimizers per
#             protein, k 7 ladder (window 4) for leaves with < 32 informative 11-mers
#   <alphabet>  the same plus one layer on every leaf: the rarest k-mer per window of 8 in that
#             alphabet, k chosen so every layer has ~28 chance-match bits per k-mer
#             (alphabets.py --budgets 28; the quarter sample's ~20 bits plus log2 of the
#             subset's ~200x more proteins, so chance hits per key stay comparable)
# The old recall tier (Dayhoff k 7, ~15 bits) is left out: at this size its query would expand
# ~200x the quarter sample's ~2.5e9 rows.
# a leaf per distinct dense 8-mer set, the dense 8-mer topology, Fitch diffs on the spanning
# tree itself. PARTITION (~130 GB, one job) splits ~53 M proteins into parts of 2 M; BUILD runs
# per part and tier (~60 GB with a layer: 7.3 GB per 262 k proteins locally); MERGE writes one
# index per tier. TARGETS/SIMULATE: 2,000 proteins, one strain each at 85/80/75/70% identity
# (BLOSUM62, depth lognormal mu 0.3 sigma 0.6); QUERY per tier with decoys (reads reversed)
# setting a 1% chance-share bar; SCORE and SUMMARY.
# Read:
#   results/bench.tsv                per tier: postings (x_base: over standard), decoy bar and
#                                    seeds, candidates, recall at depth > 1 by identity
#                                    (ev_0.01_nb: seeds at the decoy bar, neighbour rule)
#   results/index-<tier>/meta.json   postings by kind, bytes, positions, parts, build time
#   results/runs/s-<tier>/           recall.tsv (identity x depth), strategies.tsv
#   results/trace.tsv                time and peak RSS per task
# Re-submitting resumes: PARTITION is reused; BUILD and MERGE rerun when build.py changes.
source ../common.sh
[[ -f $SUBSET/members.parquet ]] || { echo "no $SUBSET/members.parquet" >&2; exit 1; }
layer() { echo "$1:--layers $2:8:rare"; }
tiers=(
    "standard:--mode standard"
    "$(layer protein protein:7)"
    "$(layer opt12 opt12:8)"
    "$(layer seb14 seb14:8)"
    "$(layer opt10 opt10:9)"
    "$(layer murphy15 murphy15:9)"
    "$(layer opt8 opt8:9)"
    "$(layer murphy10 murphy10:10)"
    "$(layer seb8 seb8:11)"
    "$(layer murphy8 murphy8:11)"
    "$(layer opt6 opt6:11)"
    "$(layer diamond11 diamond11:12)"
    "$(layer dayhoff dayhoff:13)"
    "$(layer opt4 opt4:15)"
)
TIERS=$(IFS=';'; echo "${tiers[*]}")
nextflow run "$WF/tree-index/build.nf" -c ../site.config -resume --outdir results \
    --members "$(realpath "$SUBSET/members.parquet")" --tiers "$TIERS"
