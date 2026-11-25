#!/usr/bin/env bash
#SBATCH --job-name=kfp-anchor-large
#SBATCH --output=%x.%j.log
#SBATCH --export=ALL
# The MGnify90 residue scan (workflows/pfam-anchor/lookup) on one very large Pfam family at full release size: does recall
# hold past the 58,000 rows tested on the laptop, and what do the scan and the selection of the top rows cost on x86?
#   FAMILY=5 PARQUETS='/path/results_main_<arm>/*/function-summary/pfam/*_pfam_alignments.parquet' \
#   PFAM=/path/Pfam-A.hmm ./run.sh submit
# Settings (environment):
#   FAMILY       accession number: 5 is PF00005 (62.6 M domains in the release), 72 is PF00072 (22.5 M)
#   PARQUETS     the raw-reads pipeline's alignments (seeded engines), any Pfam release: the reads are re-aligned (quote it: it is a glob)
#   PFAM         Pfam-A.hmm[.gz], the release to align both sides to
#   RELEASE      release directory or URL with mgy_proteins_pfam.parquet (default: the FTP site; about 110 GB, scanned)
#   REPS         mgy_clusters.fa.gz (default: $RELEASE's). A URL is streamed (157 GB), never stored
#   PER_LIBRARY  reads per library (200)   CHUNKS  alignment tasks (30)   EXPECT  records REPS must hold (1661142385)
#   SBATCH_ARGS  extra sbatch options for every job, a partition say. A log name needs %j, or the jobs overwrite each other's
# Steps: prepare -> align (array) -> pack -> truth -> scan -> report; ./run.sh STEP runs one here.
# Read: PF<family>/report.tsv (recall by subsample, K and SECOND), report_by_library.tsv, null.tsv (best score of reversed
#   reads, by family size), size.tsv (bytes the family needs), scan*.times.tsv (core-seconds by stage).
# Estimates, not measured: prepare 2-8 h (streaming), align 5-10 CPU-h in all, truth about 40-80 CPU-h for PF00005,
#   scan 2-6 CPU-h; 60 GB of disk in PF<family>/ and 64 GB of memory at most.
# Re-running a step skips what it has finished. Needs the repo's venv (setup.sh) and singularity, and cargo or else the Rust image, pulled once by prepare.
set -Eeuo pipefail   # -E: the ERR trap fires inside functions too
trap 'echo "failed at line $LINENO: $BASH_COMMAND" >&2' ERR
export LC_ALL=C
# a job runs a copy of this script, so submit says where it is. Not SLURM_SUBMIT_DIR: submitting from an interactive
# job, that is where the interactive job was started
HERE=${KFP_RUN_DIR:-$(cd "$(dirname "$0")" && pwd)}
cd "$HERE"
REPO=$(cd "${REPO:-../../../kmer-functional-profiler}" && pwd)
LOOKUP=$REPO/workflows/pfam-anchor/lookup
PY=${PY:-$REPO/.venv/bin/python}
: "${FAMILY:?set FAMILY, 5 for PF00005}"
RELEASE=${RELEASE:-https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release}
REPS=${REPS:-$RELEASE/mgy_clusters.fa.gz}
EXPECT=${EXPECT:-1661142385}
CHUNKS=${CHUNKS:-30}
RUNTIME=${RUNTIME:-singularity}
set_threads () { export THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-4}}; }
D=$HERE/$(printf 'PF%05d' "$FAMILY")
mkdir -p "$D"
large () { "$PY" "$LOOKUP/large.py" "$@" --out "$D"; }

# Built on the node it runs on, for that node's CPU: an AVX-512 binary would crash on an older one. With cargo, or in the
# Rust image, which the scan then runs in too, as a binary of its glibc may not run on the node. Each job builds in its
# own local temp directory (about a minute): jobs sharing a target directory on a network filesystem clobber each
# other's objects, as cargo's lock does not hold across nodes.
command -v cargo > /dev/null || [[ ! -x $HOME/.cargo/bin/cargo ]] || PATH=$HOME/.cargo/bin:$PATH
RUST_SIF=${RUST_SIF:-$HERE/rust.sif}
build () {
    T=$(mktemp -d "${TMPDIR:-/tmp}/kfp-lookup.XXXXXX")
    trap 'rm -rf "$T"' EXIT
    export RUSTFLAGS='-C target-cpu=native' CARGO_TARGET_DIR=$T/target CARGO_HOME=$T/cargo
    if command -v cargo > /dev/null; then
        cargo build --quiet --release --offline --manifest-path "$LOOKUP/rs/Cargo.toml"
        lookup () { "$CARGO_TARGET_DIR/release/lookup" "$@"; }
    else
        if [[ ! -s $RUST_SIF ]]; then
            singularity pull "$RUST_SIF.$$" "${RUST_IMAGE:-docker://rust:1.92-slim-bookworm}"
            mv "$RUST_SIF.$$" "$RUST_SIF"
        fi
        in_sif () { singularity exec -B "$REPO,$D,$T" "$RUST_SIF" "$@"; }
        # its own toolchain by path: rustup would fetch the one rust-toolchain.toml names
        # cargo's settings passed explicitly: the image's ENV CARGO_HOME (read-only in a sif) can win over the host's
        in_sif env CARGO_HOME="$CARGO_HOME" CARGO_TARGET_DIR="$CARGO_TARGET_DIR" RUSTFLAGS="$RUSTFLAGS" \
            sh -c 'PATH=$(echo /usr/local/rustup/toolchains/*/bin):$PATH; cargo build --quiet --release --offline --manifest-path "$1"' \
            sh "$LOOKUP/rs/Cargo.toml"
        lookup () { in_sif "$CARGO_TARGET_DIR/release/lookup" "$@"; }
    fi
}

[[ ${1:-} == submit ]] || set_threads   # not in submit, whose jobs inherit its environment

case ${1:-} in
submit)
    : "${PARQUETS:?set PARQUETS}" "${PFAM:?set PFAM}"
    [[ -x $PY ]] || { echo "no $PY: run $REPO/workflows/setup.sh first" >&2; exit 1; }
    [[ -s $PFAM ]] || { echo "no $PFAM" >&2; exit 1; }
    compgen -G "$PARQUETS" > /dev/null || { echo "nothing matches $PARQUETS" >&2; exit 1; }
    if [[ $RELEASE != http* ]]; then
        [[ -s $RELEASE/mgy_proteins_pfam.parquet ]] || { echo "no $RELEASE/mgy_proteins_pfam.parquet" >&2; exit 1; }
    fi
    [[ $REPS == http* || -s $REPS ]] || { echo "no $REPS" >&2; exit 1; }
    [[ -s $LOOKUP/large.py ]] || { echo "no $LOOKUP/large.py: is the branch checked out in $REPO?" >&2; exit 1; }
    command -v cargo > /dev/null || command -v singularity > /dev/null || { echo "no cargo and no singularity on PATH (the jobs build the scan)" >&2; exit 1; }
    export KFP_RUN_DIR=$HERE
    # one sbatch refused (an option the site lacks, say) cancels those already queued, or a re-submit runs them twice
    queued=()
    sb () {   # sets $id
        id=$(sbatch --parsable --chdir="$HERE" ${SBATCH_ARGS:-} "$@") || {
            (( ${#queued[@]} )) && scancel "${queued[@]}"
            echo "sbatch $* failed; cancelled ${queued[*]:-nothing}" >&2; exit 1; }
        id=${id%%;*}; queued+=("$id"); }
    sb --cpus-per-task=4 --mem=48G --time=1-00:00:00 "$HERE/run.sh" prepare; j=$id
    sb --dependency="afterok:$j" --array="0-$((CHUNKS - 1))" --cpus-per-task=4 --mem=24G --time=12:00:00 "$HERE/run.sh" align; a=$id
    sb --dependency="afterok:$a" --cpus-per-task=2 --mem=48G --time=12:00:00 "$HERE/run.sh" pack; p=$id
    sb --dependency="afterok:$p" --cpus-per-task=48 --mem=48G --time=12:00:00 "$HERE/run.sh" truth; t=$id
    # timed, but not --exclusive (some sites refuse it): node neighbours make scan*.times.tsv noisier, not the recall
    sb --dependency="afterok:$p" --cpus-per-task=16 --mem=64G --time=12:00:00 "$HERE/run.sh" scan; s=$id
    sb --dependency="afterok:$t:$s" --cpus-per-task=1 --mem=16G --time=1:00:00 "$HERE/run.sh" report; r=$id
    echo "prepare $j, align $a, pack $p, truth $t, scan $s, report $r" ;;
prepare)
    : "${PARQUETS:?set PARQUETS}" "${PFAM:?set PFAM}"
    large hmm --pfam "$PFAM" --family "$FAMILY"
    large queries --parquets "$PARQUETS" --family "$FAMILY" --per-library "${PER_LIBRARY:-200}" --runtime "$RUNTIME"
    [[ -s $D/ids.u64 ]] || large ids --release "$RELEASE" --family "$FAMILY" --threads "$THREADS" --memory 40
    if [[ ! -s $D/picked.tsv ]]; then
        build
        cat=(cat "$REPS")
        [[ $REPS == http* ]] && cat=(curl -fsSL "$REPS")
        # no retry: a resumed stream would repeat or lose records, so the count is checked instead
        "${cat[@]}" | gzip -dc | lookup pick "$D/ids.u64" "$D/chunk" "$CHUNKS" > "$D/picked.tsv.part"
        read -r seen got _ < "$D/picked.tsv.part"
        [[ $seen == "$EXPECT" ]] || { echo "$REPS gave $seen records, expected $EXPECT" >&2; exit 1; }
        (( got > 0 )) || { echo "no representative has the family" >&2; exit 1; }
        mv "$D/picked.tsv.part" "$D/picked.tsv"
    fi
    cat "$D/picked.tsv" ;;
align) large align --chunk "${CHUNK:-${SLURM_ARRAY_TASK_ID:?an array task, or set CHUNK}}" --threads "$THREADS" --memory 20 --runtime "$RUNTIME" ;;
pack) large pack ;;
truth) build; lookup truthsub "$D/bin" ;;
scan)
    build
    export FLANK=10   # only the domain and 10 residues either side, as a self-contained index holds
    lookup scan "$D/bin" "$D/scan.tsv" "${CONFIGS:-100:200,500:200,2000:200,2000:1000}" 0,2,4,6
    export THREADS=1 LIMIT=${LIMIT:-200}
    lookup scan "$D/bin" "$D/scan1.tsv" 2000:1000 0 ;;   # the default since 2026-11-25
report) large report ;;
*) sed -n '5,/^set /{/^#/p;}' "$0"; exit 1 ;;
esac
