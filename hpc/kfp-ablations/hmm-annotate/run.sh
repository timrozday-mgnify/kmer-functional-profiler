#!/usr/bin/env bash
#SBATCH --job-name=kfp-hmm-annotate
#SBATCH --output=%x.%j.log
#SBATCH --time=14-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Every MGnify protein against an HMM collection, kfpfam's prefilter then HMMER on its candidates,
# or hmmsearch alone (workflows/hmm-annotate/README.md). Inputs, as environment variables:
#   ENGINE=kfpfam-top              kfpfam's top mode, kfpfam-fast, nail (its seeds at P 0.01, as B11's arm,
#                                  then HMMER), or hmmsearch: every model, no prefilter
#   MEMBERS=true                   representatives by ENGINE, then their members by HMMER alone, on the models the
#                                  representative nearly passes. Every protein still, at a fraction of the cost
#   MEMBERS=fast                   representatives by ENGINE (hmmsearch, or kfpfam-fast), and then kfpfam's fast mode on
#                                  every member too: HMMER has the models of either. MEMBERS=top: its top mode; nail: nail
#   MIN_CLUSTER=1                  with MEMBERS=fast, top or nail, ENGINE has the representatives of clusters this size or more only;
#                                  kfpfam or nail has the whole of the others
#   SIF=<kfpfam-filter.sif>        kfpfam's image, e54c21a or later (all engines); ../images/kfpfam-filter.{tar,sif} by default
#   COLLECTION=pfam PFAM=39.0      HMM=<Pfam-A.hmm[.gz]>: the FTP site's of that release if not set. 39.0 is
#                                  what the B11 arms search ($B11_DIR/pfam/39.0/Pfam-A.hmm is the same file):
#                                  the index must be of the models the reads are aligned to
#   COLLECTION=pfam PFAM=38.2      also compares with the release's own Pfam table, made with 38.2
#                                  (COMPARE=true or false says otherwise)
#   COLLECTION=kofam               HMM=<all profiles in one file>  THRESHOLDS=<ko_list[.gz]>; labels only
#                                  (with SAMPLE > 1 domains too, for the check). From KEGG's FTP site:
#                                  cat profiles/K*.hmm > kofam.hmm
#   SAMPLE=1000                    1 in N proteins, and full hmmsearch on 1 in CHECK (10) of those
#                                  (none for ENGINE=hmmsearch, but its MEMBERS)
#   RELEASE=<mirror>               a local mirror instead of the FTP site. Worth having: even 1 in 1,000
#                                  reads most of the release's 1 TB of sequences, as samples are spread out
# First, and look at check.tsv and trace.tsv before the whole:
#   PFAM=38.2 sbatch run.sh
#   PFAM=39.0 sbatch --dependency=afterany:<that job> run.sh
#   PFAM=39.0 ENGINE=kfpfam-fast sbatch --dependency=afterany:<that job> run.sh
#   PFAM=39.0 ENGINE=hmmsearch sbatch --dependency=afterany:<that job> run.sh      # the benchmark
#   PFAM=39.0 ENGINE=hmmsearch MEMBERS=fast sbatch --dependency=afterany:<that job> run.sh     # check.members.tsv
#   PFAM=39.0 ENGINE=hmmsearch SAMPLE=1 sbatch run.sh                              # or a kfpfam mode
# kfpfam's and nail's speed against hmmsearch's on the HPC's CPUs, and what it misses: benchmark.sh
# One at a time in this directory: they share its work/ and -resume cache.
# Read results-<collection><release>-<engine>[-members[-fast|-top|-nail][-min<n>]]-<sample>/ (results-pfam39.0-hmmsearch-1000):
#   check.members.tsv  with MEMBERS, the members' recall; check.tsv is the representatives'
#   check.tsv       kfpfam's recall of the passing domains of full hmmsearch, all and by model
#   compare.tsv     the same against the release's own table (38.2)
#   summary.tsv     proteins, residues, pairs, domains passing
#   trace.tsv       time and memory per task: realtime x cpus of SEARCH x SAMPLE is the whole run's cost
#   domains/        lookup-index's input (Pfam)
#   labels/         protein_id, model (KOfam: lookup-index's LABELS)
# Cost, estimates to be replaced by the 1-in-1,000 run's trace:
#   SAMPLE=1 pfam   on the laptop top was 1.6 times hmmsearch's speed, fast 9 times (README);
#                   x86 not measured. hmmsearch would be about 113,000 core-hours at
#                   laptop speed (pfam-reannotate's 63,600 for 38.2's 30,134 models, 39.0 has 53,275)
#   check           about 10 core-hours at SAMPLE=1000 CHECK=10 (hmmsearch, laptop speed)
#   disk            domains about 64 bytes a row (laptop), so some 250 GB for Pfam over every protein
#                   at 38.2's hits, more at 39.0;
#                   work/ about 1.1 TB of FASTA, deletable when the run has finished
# Re-submitting resumes.
source ../common.sh
SIF=${SIF:-$(sif kfpfam-filter)}
collection=${COLLECTION:-pfam}
sample=${SAMPLE:-1000}
engine=${ENGINE:-kfpfam-top}
check=${CHECK:-10}
members=${MEMBERS:-false}
min_cluster=${MIN_CLUSTER:-1}
case $members in
    false | true | fast | top | nail) ;;
    *) echo "MEMBERS=$members: false, true, fast, top or nail" >&2; exit 1 ;;
esac
case $engine in
    hmmsearch | nail | kfpfam-top | kfpfam-fast) ;;
    *) echo "ENGINE=$engine: hmmsearch, nail, kfpfam-top or kfpfam-fast" >&2; exit 1 ;;
esac
filter=
case $members in
    fast | top) filter=kfpfam ;;
    nail) filter=nail ;;
esac
if [[ $filter == kfpfam && $engine == kfpfam-* && $engine != "kfpfam-$members" ]]; then
    echo "MEMBERS=$members: not ENGINE=$engine, as kfpfam has one mode to a run" >&2
    exit 1
fi
if [[ $min_cluster != 1 && -z $filter ]]; then
    echo "MIN_CLUSTER needs MEMBERS=fast, top or nail" >&2
    exit 1
fi
if [[ $sample == 1 || ( $engine == hmmsearch && $members == false ) ]]; then check=0; fi
args=(--engine "${engine%%-*}")
if [[ $engine == kfpfam-* ]]; then args+=(--kfpfam_mode "${engine#*-}"); fi
name=$collection
case $collection in
    pfam)
        release=${PFAM:-39.0}
        name=pfam$release
        if [[ $members == false ]]; then
            args+=(--compare_release "${COMPARE:-$([[ $release == 38.2 ]] && echo true || echo false)}")
        fi
        args+=(--hmm "${HMM:-https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam$release/Pfam-A.hmm.gz}")
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
targets=all
if [[ $members != false ]]; then
    targets=reps engine=$engine-members
    # a sample's members are 1 in $sample of every range of ids: wider ones, so that a task is not mostly start-up
    args+=(--members true --member_width $((20000000 * (sample < 100 ? sample : 100))))
fi
if [[ -n $filter ]]; then
    engine+=-$members
    args+=(--member_filter "$filter" --min_cluster "$min_cluster")
    if [[ $filter == kfpfam && $engine != kfpfam-* ]]; then args+=(--kfpfam_mode "$members"); fi
    if [[ $min_cluster != 1 ]]; then engine+=-min$min_cluster; fi
fi
nextflow run "$WF/hmm-annotate" -c ../site.config \
    ${RELEASE:+--release "$RELEASE"} --kfpfam_sif "$(realpath "$SIF")" "${args[@]}" \
    --targets "$targets" --sample "$sample" --check_sample "$check" \
    --outdir "results-$name-$engine-$sample" -resume
