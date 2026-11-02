#!/usr/bin/env bash
#SBATCH --job-name=kfp-tree-index-slice
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Tree index on a family-complete slice (docs/tree-index-plan.md, Progress log: precision and
# read classification), with workflows/tree-index/build.nf. A cluster sample (the quarter
# sample, the 1-in-100 subset) almost never holds a protein's neighbouring units, so read
# precision by unit cannot be measured on it; the slice keeps every unit of chosen families.
#   FAMILIES  40 Pfam families of 20 k-150 k annotated proteins, at random (slice.py), from
#             the release's Pfam table (release-tiers' results-standard/pfam.parquet)
#   TAKE      per members bucket (256): every member of each MGnify90 cluster holding one of
#             their proteins; SLICE concatenates (~2-6 M proteins, estimate)
#   then the tiers (slim, standard, large) built, queried and aligned as tree-index-build's,
#   with recall.py reads: read sensitivity and precision by leaf, group, unit and family.
# Read:
#   results/families.tsv, families.json  the families chosen and their proteins
#   results/members.parquet              the slice
#   results/bench.tsv                    per tier: recall by identity, candidate unit
#                                        precision, read sensitivity and precision (unit, Pfam)
#   results/runs/s-<tier>/reads.tsv      per identity bin: reads, called, ties, precision by
#                                        leaf / group / unit / pfam, sensitivity
#   results/runs/s-<tier>/reads.json     the same overall, candidate precision, read bar
#   results/trace.tsv                    time and peak RSS per task
# Re-submitting resumes.
source ../common.sh
STD=../release-tiers/results-standard
compgen -G "$STD/members.bucket*.parquet" >/dev/null || { echo "no $STD/members.bucket*.parquet: run release-tiers first" >&2; exit 1; }
[[ -f $STD/pfam.parquet ]] || { echo "no $STD/pfam.parquet (release-tiers with --pfam)" >&2; exit 1; }
nextflow run "$WF/tree-index/build.nf" -c ../site.config -resume --outdir results \
    --slice_store "$(realpath "$STD")/members.bucket*.parquet" \
    --slice_pfam "$(realpath "$STD/pfam.parquet")" \
    --tiers 'slim:--mode slim;standard:--mode standard;large:--mode large'
