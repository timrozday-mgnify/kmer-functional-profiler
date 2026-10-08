#!/usr/bin/env bash
#SBATCH --job-name=kfp-release-tiers
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# The three release tiers over the whole MGnify Proteins release (plan, decided 2026-10-06),
# each with mgnify-subset's partitioned build, into results-<tier>/index:
#   slim      t_base 0.01, singletons 0.001, n_min 8               ~91 GB today's format
#   standard  tbase02_all16: t_base 0.02, n_min 16, oversample 1  ~196 GB
#   large     tbase05_all: t_base 0.05, n_min 32, oversample 1    ~282 GB
# (sizes from the 1-in-100 builds x 100, +-20-30%). One after another in this directory, so
# -resume reuses the release extraction (MEMBERSHIP, EXTRACT, MERGE: ~1.3 TB read, ~1 TB of
# cluster buckets linked under each results-<tier>/) and only the build runs again.
# RELEASE=<local mirror> sbatch run.sh reads a mirror instead of the FTP site.
# Each index gets the aai survival model step34q attached (plan, decided 2026-10-08).
for tier in "slim:--t-base 0.01 --t-base-singleton 0.001 --n-min 8:128" \
            "standard:--t-base 0.02 --n-min 16 --oversample 1:256" \
            "large:--t-base 0.05 --n-min 32 --oversample 1:256"; do
    IFS=: read -r name args concat_gb <<< "$tier"
    # a config file, not --index_args: Nextflow reads a value starting with '--' as an option.
    # CONCAT holds tier 2 (sparse at 77 GB: ~31 GB) and BLOOM the candidate filter (~14 GB
    # sparse), both larger at higher floors: give them memory up front, not through retries.
    # PRESENCE, GROUPS, PACK_RANGE (32 GB) and UNITS (128 GB) retried on the large tier.
    cat > "$name.config" <<CFG
params.index_args = "--k 11 $args"
process {
    withName: 'CONCAT' { memory = { ${concat_gb}.GB * task.attempt } }
    withName: 'BLOOM' { memory = { 64.GB * task.attempt } }
    withName: 'PRESENCE|GROUPS|PACK_RANGE' { memory = { 64.GB * task.attempt } }
    withName: 'UNITS' { memory = { 384.GB * task.attempt } }
}
CFG
    nextflow run "$WF/mgnify-subset" -c ../site.config -c "$name.config" ${RELEASE:+--release "$RELEASE"} \
        --biome root --sample 1 --shards 256 --buckets 256 --ranges 256 --index false --build true \
        --publish_mode link --outdir "results-$name" -resume
    "$REPO/.venv/bin/kmer-functional-profiler" aai-model "results-$name/index" \
        "$WF/fmh-benchmark/ablations/aai_models/step34q.json"
done
