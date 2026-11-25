#!/usr/bin/env bash
#SBATCH --job-name=kfp-hmm-annotate
#SBATCH --output=%x.%j.log
#SBATCH --time=14-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Every MGnify protein against an HMM collection, kfpfam's prefilter then HMMER on its candidates
# (workflows/hmm-annotate/README.md). Inputs, as environment variables:
#   SIF=<kfpfam-filter.sif>        kfpfam's image, 86f2a63 or later; ../images/kfpfam-filter.{tar,sif} by default
#   COLLECTION=pfam                HMM=<Pfam-A.hmm[.gz]>: 38.2 from the FTP site if not set. The
#                                  index must be of the models the reads are aligned to
#   COLLECTION=kofam               HMM=<all profiles in one file>  THRESHOLDS=<ko_list[.gz]>; labels only
#                                  (with SAMPLE > 1 domains too, for the check). From KEGG's FTP site:
#                                  cat profiles/K*.hmm > kofam.hmm
#   SAMPLE=1000                    1 in N proteins, and full hmmsearch on 1 in CHECK (10) of those
#   RELEASE=<mirror>               a local mirror instead of the FTP site. Worth having: even 1 in 1,000
#                                  reads most of the release's 1 TB of sequences, as samples are spread out
# First, and look at check.tsv and trace.tsv before the whole:
#   COLLECTION=pfam sbatch run.sh
#   COLLECTION=pfam SAMPLE=1 sbatch run.sh
# Read results-<collection>-<sample>/:
#   check.tsv       kfpfam's recall of the passing domains, all and by model
#   summary.tsv     proteins, residues, pairs, domains passing
#   trace.tsv       time and memory per task: realtime x cpus of SEARCH x SAMPLE is the whole run's cost
#   domains/        lookup-index's input (Pfam)
#   labels/         protein_id, model (KOfam: lookup-index's LABELS)
# Cost, estimates to be replaced by the 1-in-1,000 run's trace:
#   SAMPLE=1 pfam   not measured on proteins. hmmsearch alone would be about 63,600 core-hours at
#                   laptop speed (pfam-reannotate's README); kfpfam on reads' ORFs was well under nail
#                   there but about 13x slower per read on Codon
#   check           about 6 core-hours at SAMPLE=1000 CHECK=10 (hmmsearch, laptop speed)
#   disk            domains about 64 bytes a row (laptop), so some 250 GB for Pfam over every protein;
#                   work/ about 1.1 TB of FASTA, deletable when the run has finished
# Re-submitting resumes.
source ../common.sh
SIF=${SIF:-$(sif kfpfam-filter)}
collection=${COLLECTION:-pfam}
sample=${SAMPLE:-1000}
check=$([[ $sample == 1 ]] && echo 0 || echo "${CHECK:-10}")
args=()
case $collection in
    pfam)
        args+=(--hmm "${HMM:-https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/Pfam-A.hmm.gz}")
        ;;
    kofam)
        : "${HMM:?KOfam profiles in one file}" "${THRESHOLDS:?KOfam ko_list}"
        args+=(--hmm "$HMM" --thresholds "$THRESHOLDS")
        if [[ $sample == 1 ]]; then args+=(--domains false); fi
        ;;
    *)
        : "${HMM:?the collection as one HMM file}"
        args+=(--hmm "$HMM" ${THRESHOLDS:+--thresholds "$THRESHOLDS"})
        ;;
esac
nextflow run "$WF/hmm-annotate" -c ../site.config \
    ${RELEASE:+--release "$RELEASE"} --kfpfam_sif "$(realpath "$SIF")" "${args[@]}" \
    --targets all --sample "$sample" --check_sample "$check" \
    --outdir "results-$collection-$sample" -resume
