#!/usr/bin/env bash
#SBATCH --job-name=kfp-tree-index-wide
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Tree index, tier 1 as fewer representatives sensitive to wider AAI (docs/tree-index-plan.md,
# Progress log: layer-only tiers, wide clusters for key budgets), on the 1-in-100 subset with
# workflows/tree-index/build.nf, confirming the mini-slice grid at scale. Layer only (no
# per-leaf base: --t-base 0 --window 0 --ladder= --n-1 0); keys on each tree group's
# composite representative, budgets spent per wide cluster of groups (representatives sharing
# >= 30% of 5-mers, ~79% identity: --layer-units kmer:5:0.3), window 8:
#   wide-k4-3x   layer k 4, ~3x today's keys (sqrt:9.8)   mini-slice: 0.96 / 0.93 / 0.95 at 80/75/70%
#   wide-k4-4x   layer k 4, ~4x (sqrt:13)                 0.96 / 0.96 / 0.97
#   wide-k5-6x   layer k 5, ~6x (sqrt:19.5)               0.96 / 0.94 / 0.94
#   wide-k5-8x   layer k 5, ~8x (sqrt:26)                 0.99 / 0.95 / 0.97
#   large-k7     the full-leaf large tier at layer k 7: the control (0.89 / 0.79 / 0.62 here)
# --chance-scale 108 (MGnify proteins / the subset's): decoy counts scaled so the bars are the
# release's. Short layer k aligns many chance pairs (k 4: ~1,500 per read on the mini-slice,
# x the database): a tier whose QUERY fails twice is skipped, so the rest still score.
# Read:
#   results/bench.tsv          per tier: recall by identity (verify_0.01_nb), release sizes
#                              (release_lookup_gb with positions per leaf, as built;
#                              release_align_gb), pairs per read now and at release
#   results/runs/q-<tier>/     query.json, verify.json (bars, pairs aligned), tier2/
#   results/runs/s-<tier>/     strains.parquet, reads.tsv, detection.tsv (tier 1 and tier2/)
#   results/trace.tsv          time and peak RSS per task (QUERY for k 4 especially)
# Re-submitting resumes.
source ../common.sh
[[ -f $SUBSET/members.parquet ]] || { echo "no $SUBSET/members.parquet" >&2; exit 1; }
LO="--t-base 0 --window 0 --ladder= --n-1 0 --layer-units kmer:5:0.3"
tiers=(
    "wide-k4-3x:--mode standard $LO --layers protein:4:8:rare --layer-keys sqrt:9.8"
    "wide-k4-4x:--mode standard $LO --layers protein:4:8:rare --layer-keys sqrt:13"
    "wide-k5-6x:--mode standard $LO --layers protein:5:8:rare --layer-keys sqrt:19.5"
    "wide-k5-8x:--mode standard $LO --layers protein:5:8:rare --layer-keys sqrt:26"
    "large-k7:--mode large --layers protein:7:8:rare"
)
TIERS=$(IFS=';'; echo "${tiers[*]}")
cat > build.config <<'CFG'
process {
    withName: '.*BUILD' { memory = { 128.GB * task.attempt } }
    withName: '.*QUERY' { errorStrategy = { task.attempt < 2 ? 'retry' : 'ignore' } }
}
CFG
nextflow run "$WF/tree-index/build.nf" -c ../site.config -c build.config -resume --outdir results \
    --members "$(realpath "$SUBSET/members.parquet")" --pfam "$(realpath "$SUBSET/pfam.parquet")" \
    --tiers "$TIERS" --part_proteins 1000000 --chance_scale 108
