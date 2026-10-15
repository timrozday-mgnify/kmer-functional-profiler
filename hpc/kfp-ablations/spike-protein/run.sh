#!/usr/bin/env bash
#SBATCH --job-name=kfp-spike-protein
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Spike-in benchmark, amino-acid indexes (plan, phase 13: steps 3-6 and the crude tool's gate;
# step 6's aai check with truth). 2000 families (3 units of 3 members) from the 1-in-100
# subset's representatives; 10 samples, each spiked into the first 10 M pairs of ERR7738575
# (Hadza gut); 10% of units with a withheld sibling. One arm per release tier, step34q
# attached; DIAMOND follow-up. Read per tier results-<tier>/summary.tsv:
#   read_precision >= 0.95 and sibling_caught >= 0.8 (gate; _parent_* says why it is not),
#   aln_gap_*/aln_cover (step 6's rule), sim_* (aai against truth), f1/l1;
#   trace.tsv: PILEUP:DIAMOND against PILEUP:QUERY realtime (step 3's rule: <= ~5x),
#   <sample>/stats.json: hits_per_pair (step 0). Needs nothing from other runs.
source ../common.sh
fetch_pairs bg 10000000 ERR7738575
[[ -s seeds.faa ]] || "$REPO/.venv/bin/python" - "$SUBSET/members.parquet" seeds.faa <<'PY'
import sys
import polars as pl
reps = (
    pl.scan_parquet(sys.argv[1])
    .filter(pl.col("protein_id").cast(pl.String) == pl.col("cluster_rep").cast(pl.String),
            pl.col("full_length"), pl.col("sequence").str.len_bytes().is_between(200, 600))
    .select("protein_id", "sequence").head(2000).collect()
)
with open(sys.argv[2], "w") as f:
    f.writelines(f">{i}\n{s}\n" for i, s in reps.iter_rows())
PY
cat > spike.config <<CFG
process {
    withName: 'SIMULATE' { memory = 16.GB; time = 8.h }
    withName: 'PILEUP:QUERY' { cpus = 8; memory = { 32.GB * task.attempt } }
    withName: 'PILEUP:DIAMOND' { cpus = 16; memory = { 64.GB * task.attempt }; time = { 48.h * task.attempt } }
    withName: 'PILEUP:CALL' { memory = { 64.GB * task.attempt } }
    withName: 'PILEUP:MAP_CONSENSUS' { cpus = 16; memory = { 32.GB * task.attempt }; time = { 24.h * task.attempt } }
    withName: 'PILEUP:(HANDOFF|BAM|FREEBAYES|CONSEQUENCES)' { memory = { 16.GB * task.attempt }; time = { 48.h * task.attempt } }
    withName: 'PILEUP:(FLORIA|INSTRAIN|CLIQUESNV)' { cpus = 16; memory = { 64.GB * task.attempt }; time = { 48.h * task.attempt } }
}
CFG
# one after another in this directory, so -resume reuses the families, siblings and samples
for tier in "slim:--t-base 0.01 --t-base-singleton 0.001 --n-min 8" \
            "standard:--t-base 0.02 --n-min 16 --oversample 1" \
            "large:--t-base 0.05 --n-min 32 --oversample 1"; do
    IFS=: read -r name args <<< "$tier"
    nextflow run "$WF/spike-benchmark" -c ../site.config -c spike.config -resume \
        --seed_fasta seeds.faa --families 2000 --seeds 10 \
        --siblings_args '--siblings 0.1' --background "$PWD/bg_1.fastq.gz,$PWD/bg_2.fastq.gz" \
        --index_args "--k 11 $args" --aai_model "$WF/fmh-benchmark/ablations/aai_models/step34q.json" \
        --query_stats true --outdir "results-$name"
done
