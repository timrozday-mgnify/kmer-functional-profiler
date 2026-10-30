#!/usr/bin/env bash
#SBATCH --job-name=kfp-tree-index-build
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Tree index on the 1-in-100 subset (docs/tree-index-plan.md, Progress log: the standard index
# restored), with workflows/tree-index/build.nf:
#   standard  build.py --mode standard: t_base 0.01, no floors, ~4 spread 11-mer minimizers per
#             protein, k 7 ladder (window 4) for leaves with < 32 informative 11-mers
#   recall    build.py --mode recall: the same plus the rarest Dayhoff 7-mer per window of 8 on
#             every leaf (high recall at 70-85% identity; not yet tested for precision or F1)
# a leaf per distinct dense 8-mer set, the dense 8-mer topology, Fitch diffs on the spanning
# tree itself. The recall tier's BUILD needs more memory per part (the Dayhoff frequency table
# and per-k-mer hashing): on phase 0's quarter sample 9.4 GB per 262 k proteins. PARTITION (~130 GB, one job)
# splits ~53 M proteins into parts of 2 M; BUILD runs per part and tier (~26 GB each); MERGE
# writes one index per tier.
# Read:
#   results/index-<tier>/meta.json   postings by kind, bytes, positions, parts, build time
#   results/partition.json           components, the largest, part sizes
#   results/trace.tsv                time and peak RSS per task
# Re-submitting resumes: PARTITION is reused; BUILD and MERGE rerun when build.py changes.
source ../common.sh
[[ -f $SUBSET/members.parquet ]] || { echo "no $SUBSET/members.parquet" >&2; exit 1; }
nextflow run "$WF/tree-index/build.nf" -c ../site.config -resume --outdir results \
    --members "$(realpath "$SUBSET/members.parquet")"
