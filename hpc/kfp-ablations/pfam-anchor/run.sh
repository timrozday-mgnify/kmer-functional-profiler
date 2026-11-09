#!/usr/bin/env bash
#SBATCH --job-name=kfp-pfam-anchor
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# The Pfam-anchored route and its benchmark (docs/pfam-anchor-plan.md, phase 5) on the
# tree-index slice (40 families, every unit of their clusters), with workflows/pfam-anchor.
#   sbatch run.sh                          hmmsearch hits only
#   KFPFAM=<kfpfam checkout> sbatch run.sh   also kfpfam's weak-read recruitment (needs its
#                                          data/prep/models.tsv)
#   HMM=<Pfam-A.hmm[.gz]> sbatch run.sh    a local Pfam-A file instead of downloading it
# Read:
#   results/report.tsv            hit recall by identity bin
#   results/oracle.tsv, alignment.tsv   the oracle on every read hit
#   results/phase3_*.tsv          the anchored route against DIAMOND
#   results/mapall_*.tsv, truthscore_*.tsv   mapping, gather and pileups against the truth
#   results/trace.tsv             time and peak RSS per task (the largest families' cost)
#   results/data/                 every intermediate (the scripts' layout)
# Re-submitting resumes.
source ../common.sh
SLICE=../tree-index-slice/results
[[ -f $SLICE/members.parquet && -f $SLICE/families.tsv ]] || { echo "no $SLICE/members.parquet: run tree-index-slice first" >&2; exit 1; }
if [[ -z ${HMM:-} ]]; then
    HMM=Pfam-A.hmm.gz
    [[ -s $HMM ]] || curl -fsSL --retry 5 -o "$HMM.tmp" \
        https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/Pfam-A.hmm.gz && mv "$HMM.tmp" "$HMM"
fi
nextflow run "$WF/pfam-anchor" -c ../site.config -profile singularity -resume --outdir results \
    --members "$(realpath "$SLICE/members.parquet")" --families "$(realpath "$SLICE/families.tsv")" \
    --pfam "$(realpath "$HMM")" ${KFPFAM:+--kfpfam "$(realpath "$KFPFAM")"} \
    --recruit "$([[ -n ${KFPFAM:-} ]] && echo true || echo false)"
