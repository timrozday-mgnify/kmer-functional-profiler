#!/usr/bin/env bash
#SBATCH --job-name=kfp-lookup-profile
#SBATCH --output=%x.%A_%a.log
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --export=ALL
# Profiles Pfam alignments that exist already (*_pfam_alignments.parquet, the raw-reads pipeline's nail and
# kfpfam engines) with the lookup bundle, an array task each: what its PFAM_LOOKUP does, without re-running it.
#   HMM=$B11_DIR/pfam/39.0/Pfam-A.hmm \
#   PARQUETS="$B11_DIR/pfam39.0/results_main_kfpfam_a/*/function-summary/pfam/*_pfam_alignments.parquet" ./run.sh
#   HMM       the unzipped Pfam-A.hmm the reads were aligned to. It must be the bundle's, by SHA-256: that is the only
#             check of the release on alignments older than the pfam_version column
#   PARQUETS  a glob, quoted, or a list: one arm's alignments, a file per sample
#   BUNDLE=../lookup-index/results/bundle    OUT=results
#   SIF       ../images/kfp-lookup.{tar,sif} by default
#   FLOOR=0.4 K=2000 SECOND=1000 MARGIN=40   bundle.py's settings, as the pipeline's PFAM_LOOKUP (MARGIN=0: none)
#   MEM=128G  if the largest family does not fit: the index holds one family at a time (PF00005 30-40 GB, estimated)
#   DRY_RUN=1 lists the samples and submits nothing
# Read, in OUT:
#   <sample>_pfam_lookup.reads.parquet   every alignment's call: protein (its cluster's representative), similarity, labels
#   <sample>_pfam_lookup.ko.tsv          reads per label
#   <sample>_pfam_lookup.summary.json    alignments, reads called, under the floor
# Re-running submits only the samples with no summary.json.
source ../common.sh
OUT=${OUT:-results}
if [[ -z ${SLURM_ARRAY_TASK_ID:-} ]]; then
    : "${HMM:?the unzipped Pfam-A.hmm the reads were aligned to}" "${PARQUETS:?the alignments, a quoted glob}"
    BUNDLE=$(realpath "${BUNDLE:-../lookup-index/results/bundle}") HMM=$(realpath "$HMM")
    [[ -s $BUNDLE/manifest.json && -s $HMM ]] || { echo "no $BUNDLE/manifest.json or $HMM" >&2; exit 1; }
    SIF=${SIF:-$(sif kfp-lookup)} && SIF=$(realpath "$SIF")
    mkdir -p "$OUT"
    for p in $PARQUETS; do
        [[ -s $p ]] || { echo "no $p" >&2; exit 1; }
        s=$(basename "$p" _pfam_alignments.parquet)
        [[ -s $OUT/${s}_pfam_lookup.summary.json ]] || printf '%s\t%s\n' "$s" "$(realpath "$p")"
    done > "$OUT/samples.tsv"
    twice=$(cut -f1 "$OUT/samples.tsv" | sort | uniq -d)
    [[ -z $twice ]] || { echo "in two files, give one arm's alignments only:" $twice >&2; exit 1; }
    n=$(( $(wc -l < "$OUT/samples.tsv") ))
    cat "$OUT/samples.tsv"
    (( n > 0 )) || { echo "nothing to do"; exit 0; }
    [[ ${DRY_RUN:-0} == 1 ]] && exit 0
    export BUNDLE HMM SIF OUT
    sbatch --array="1-$n" ${MEM:+--mem="$MEM"} run.sh
    exit
fi
IFS=$'\t' read -r sample parquet < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$OUT/samples.tsv")
out=$(realpath "$OUT")
singularity exec --bind "$BUNDLE" --bind "$(dirname "$HMM")" --bind "$(dirname "$parquet")" --bind "$out" "$SIF" \
    python3 /opt/lookup/bundle.py profile --bundle "$BUNDLE" --hmm "$HMM" --parquet "$parquet" \
    --out "$out/${sample}_pfam_lookup" --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --floor "${FLOOR:-0.4}" --k "${K:-2000}" --second "${SECOND:-1000}" --margin "${MARGIN:-40}"
