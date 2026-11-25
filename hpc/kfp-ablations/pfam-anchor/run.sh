#!/usr/bin/env bash
#SBATCH --job-name=kfp-pfam-anchor
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# The Pfam-anchored route and its benchmark (docs/pfam-anchor-plan.md, phases 5 and 6) on the
# tree-index slice (40 families, every unit of their clusters), with workflows/pfam-anchor.
# Synthetic reads on a depth x AAI grid (strains at 100-70% identity, held-out relatives).
#   sbatch run.sh                          hmmsearch hits only
#   KFPFAM=<kfpfam checkout> sbatch run.sh   also kfpfam's weak-read recruitment (needs its
#                                          data/prep/models.tsv)
#   HMM=<Pfam-A.hmm[.gz]> sbatch run.sh    a local Pfam-A file instead of downloading it
# The same spikes in real metagenomes: pfam-anchor-{gut,hadza,soil}/run.sh (BG=<run>), after
# this run's annotation (they reuse results/data's hits); BG_PAIRS background pairs (2 M).
# Read:
#   results/bench_recall.tsv      per depth x AAI cell: recall over pairs and strains (Wilson
#                                 intervals), recall against DIAMOND-resolvable pairs, ceiling
#   results/bench_misses.tsv      the stage each true pair was lost at, per cell, with examples
#   results/bench_units.tsv       unit precision (false units; background-supported units)
#   results/bench_background.tsv  real background reads against DIAMOND (with BG)
#   results/bench_cost.tsv        per stage CPU-s, wall-s, peak RSS, core-ms per pair, per Mbp;
#                                 index sizes
#   results/report.tsv            hit recall by identity bin
#   results/phase3_*.tsv, mapall_*.tsv, truthscore_*.tsv   the earlier phases' tables
#   results/trace.tsv             time and peak RSS per task
#   results/data/                 every intermediate (the scripts' layout)
# Re-submitting resumes.
source ../common.sh
SLICE=../tree-index-slice/results
[[ -f $SLICE/members.parquet && -f $SLICE/families.tsv ]] || { echo "no $SLICE/members.parquet: run tree-index-slice first" >&2; exit 1; }
if [[ -z ${HMM:-} ]]; then
    HMM=../pfam-anchor/Pfam-A.hmm.gz  # one download, shared with the BG runs
    [[ -s $HMM ]] || curl -fsSL --retry 5 -o "$HMM.tmp" \
        https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/Pfam-A.hmm.gz && mv "$HMM.tmp" "$HMM"
fi
bg=()
if [[ -n ${BG:-} ]]; then
    [[ -f ../pfam-anchor/results/data/hits.parquet ]] || { echo "run ../pfam-anchor first (its annotation is reused)" >&2; exit 1; }
    fetch_pairs bg "${BG_PAIRS:-2000000}" "$BG"
    bg=(--background_1 "$(realpath bg_1.fastq.gz)" --background_2 "$(realpath bg_2.fastq.gz)"
        --background_pairs 0 --annotated "$(realpath ../pfam-anchor/results/data)")
fi
nextflow run "$WF/pfam-anchor" -c ../site.config -profile singularity -resume --outdir results \
    --members "$(realpath "$SLICE/members.parquet")" --families "$(realpath "$SLICE/families.tsv")" \
    --pfam "$(realpath "$HMM")" ${KFPFAM:+--kfpfam "$(realpath "$KFPFAM")"} \
    --recruit "$([[ -n ${KFPFAM:-} ]] && echo true || echo false)" ${bg[@]+"${bg[@]}"}
