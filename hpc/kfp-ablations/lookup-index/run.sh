#!/usr/bin/env bash
#SBATCH --job-name=kfp-lookup-index
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# The MGnify90 lookup bundle of the whole release (workflows/lookup-index/README.md), from hmm-annotate's
# SAMPLE=1 runs of Pfam 39.0 and KOfam. Inputs, as environment variables:
#   SIF=<kfp-lookup.sif>         workflows/pfam-anchor/lookup/Dockerfile's image; ../images/kfp-lookup.{tar,sif} by default
#   HMM=<Pfam-A.hmm[.gz]>        the file hmm-annotate searched, and the one the reads are aligned to:
#                                $B11_DIR/pfam/39.0/Pfam-A.hmm, or the FTP site's 39.0 (it is checked)
#   PFAM=../hmm-annotate/results-pfam39.0-1
#   KOFAM=../hmm-annotate/results-kofam-1   its labels; KOFAM= (empty) for a bundle with none
#   SAMPLE_PARQUET=<*_pfam_alignments.parquet>   optional: profiled with the bundle at the end
#   RELEASE=<mirror>             a local mirror instead of the FTP site
#   PFAM_VERSION=39.0 MGNIFY_VERSION=2026_07     recorded in the bundle
# Read:
#   results/bundle/           the index, for the raw-reads pipeline's --pfam_lookup_bundle
#   results/bundle/manifest.json   rows, member rows, rows_label_dropped
#   results/check.*           verify, and the sample's profile
#   results/trace.tsv         time and peak memory per task
# Cost, estimates:
#   rows      about 64 bytes a domain row in SPLIT's files, about 250 GB, in some 90,000 files in work/
#   bundle    members about 53 GB (21.1 bytes x 2.5 B rows, REPORT.md), representatives about 100 bytes a row
#   memory    BUILD of the group that holds PF00005 alone 30-40 GB; asks for 64, 128 on the retry
# Re-submitting resumes.
source ../common.sh
: "${HMM:?the Pfam-A.hmm that was annotated with}"
SIF=${SIF:-$(sif kfp-lookup)}
PFAM=${PFAM:-../hmm-annotate/results-pfam39.0-1}
KOFAM=${KOFAM-../hmm-annotate/results-kofam-1}
nextflow run "$WF/lookup-index" -c ../site.config ${RELEASE:+--release "$RELEASE"} \
    --lookup_sif "$(realpath "$SIF")" --hmm "$(realpath "$HMM")" --pfam "$(realpath "$PFAM")" \
    ${KOFAM:+--labels "ko=$(realpath "$KOFAM")"} \
    ${SAMPLE_PARQUET:+--check_parquet "$(realpath "$SAMPLE_PARQUET")"} \
    --pfam_version "${PFAM_VERSION:-39.0}" --mgnify_version "${MGNIFY_VERSION:-2026_07}" \
    --outdir results -resume
