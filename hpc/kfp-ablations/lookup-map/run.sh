#!/usr/bin/env bash
#SBATCH --job-name=kfp-lookup-map
#SBATCH --output=%x.%j.log
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --export=ALL
# Every read of a sample against the MGnify proteins its lookup called, by DIAMOND (docs/pfam-anchor-plan.md,
# "Mapping every read to the lookup's proteins", step 1): the reference's size, how identical the reads are to
# it, how often a called pair's mates land on the lookup's own protein, and which seed shapes would have found
# each alignment. These choose the mapper's seeds; DIAMOND's time here is the native baseline.
#   READS=../lookup-profile/results/ERR10889147_pfam_lookup.reads.parquet RUN=ERR10889147 sbatch run.sh
#   READS     a sample's lookup calls (lookup-profile's, or the raw-reads pipeline's PFAM_LOOKUP)
#   RUN       the sample's ENA run: its first PAIRS read pairs are the queries (the lookup's reads, by name)
#   PAIRS=2000000  STORE='../release-tiers/results-standard/members.bucket*.parquet'  FLOOR=0.4
#   DIAMOND_ARGS='--sensitive -F 15 --top 10'   (the pileup workflow's)
# Read, in results/<RUN>/:
#   ref.json             proteins called, written (fewer: not in the store), residues
#   map.summary.json     reference size, mates and pairs aligned, pairs the lookup called, mates of called pairs
#                        whose best hit is the lookup's protein (the mate-first share), mates by identity, seeded
#                        share per shape
#   map.identity.tsv     the same per identity bin, called pairs apart
#   diamond.time         DIAMOND's /usr/bin/time -v, both mates (CPU-s, peak RSS)
source ../common.sh
: "${READS:?a sample's *_pfam_lookup.reads.parquet}" "${RUN:?the sample's ENA run}"
STORE=${STORE:-../release-tiers/results-standard/members.bucket*.parquet}
compgen -G "$STORE" >/dev/null || { echo "no $STORE" >&2; exit 1; }
PY=$REPO/.venv/bin/python MAPREF=$WF/pfam-anchor/lookup/mapref.py
out=results/$RUN && mkdir -p "$out"
fetch_pairs "$out/reads" "${PAIRS:-2000000}" "$RUN"
[[ -s $out/ref.faa ]] || $PY "$MAPREF" ref --reads "$READS" --store $STORE --floor "${FLOOR:-0.4}" \
    --out "$out/ref.faa" > "$out/ref.json"
cat "$out/ref.json"
sif=$(realpath ../images)/diamond-2.2.8.sif
[[ -s $sif ]] || singularity pull "$sif" docker://quay.io/biocontainers/diamond:2.2.8--he361c42_0
d() { singularity exec --bind "$(realpath "$out")" "$sif" diamond "$@"; }
t=${SLURM_CPUS_PER_TASK:-16}
d makedb --in "$out/ref.faa" --db "$out/ref" --threads "$t"
: > "$out/diamond.time"
for m in 1 2; do
    # shellcheck disable=SC2086  # DIAMOND_ARGS is several options
    /usr/bin/time -v -a -o "$out/diamond.time" singularity exec --bind "$(realpath "$out")" "$sif" \
        diamond blastx --db "$out/ref.dmnd" --query "$out/reads_$m.fastq.gz" --threads "$t" \
        ${DIAMOND_ARGS:---sensitive -F 15 --top 10} \
        --outfmt 6 qseqid sseqid pident length bitscore btop --out "$out/r$m.tsv"
done
$PY "$MAPREF" measure --reads "$READS" --ref "$out/ref.faa" --diamond "$out/r1.tsv" "$out/r2.tsv" \
    --floor "${FLOOR:-0.4}" --out "$out/map"
