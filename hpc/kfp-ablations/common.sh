# Sourced by each <run>/run.sh, which must be submitted from its own directory:
#   cd <runs>/kfp-ablations/<run> && sbatch run.sh
set -euo pipefail
REPO=../../../kmer-functional-profiler  # <runs>/../kmer-functional-profiler
WF=$REPO/workflows
[[ -f run.sh && -f ../site.config ]] || { echo "submit from the run's directory: cd <run> && sbatch run.sh" >&2; exit 1; }
[[ -x $REPO/.venv/bin/kmer-functional-profiler ]] || { echo "no venv in $REPO: run $WF/setup.sh first" >&2; exit 1; }
command -v nextflow >/dev/null || { echo "nextflow not on PATH" >&2; exit 1; }
export NXF_OPTS='-Xms1g -Xmx4g'
SUBSET=../../mgnify-proteins-subset/results/lin100-dense  # members.parquet, pfam.parquet at 1 in 100

# fmh <ablation config name, or '' for none> [nextflow options]
fmh() {
    local which=$1
    shift
    local cfg=()
    [[ -n $which ]] && cfg=(-c "$WF/fmh-benchmark/ablations/$which.config")
    nextflow run "$WF/fmh-benchmark" -c ../site.config ${cfg[@]+"${cfg[@]}"} --outdir results -resume "$@"
}

# fetch_reads RUN: an ENA run's two FASTQ files into ../reads/ (md5-checked, once, shared by
# every run directory); prints nothing. fetch_pairs OUT PAIRS RUN...: the first PAIRS read
# pairs of those runs, in order, as OUT_1.fastq.gz and OUT_2.fastq.gz (once).
fetch_reads() {
    local run=$1 urls md5s
    [[ -f ../reads/$run.ok ]] && return
    mkdir -p ../reads
    read -r _ urls md5s < <(curl -fsS --retry 5 \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$run&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" | tail -n 1)
    (cd ../reads && for url in ${urls//;/ }; do curl -fsSL --retry 5 -C - -O "https://$url"; done &&
        paste -d ' ' <(tr ';' '\n' <<< "$md5s") <(tr ';' '\n' <<< "$urls" | xargs -n 1 basename | sed 's/^/ /') | md5sum -c -)
    touch ../reads/$run.ok
}
fetch_pairs() {
    local out=$1 pairs=$2 run m
    shift 2
    [[ -s ${out}_1.fastq.gz && -s ${out}_2.fastq.gz ]] && return
    for run in "$@"; do fetch_reads "$run"; done
    for m in 1 2; do
        # head stops reading early, so the pipe's earlier stages die of SIGPIPE: no pipefail
        (set +o pipefail; for run in "$@"; do gzip -dc "../reads/${run}_$m.fastq.gz"; done |
            head -n $((4 * pairs)) | gzip -1 > "${out}_$m.fastq.gz.tmp")
        mv "${out}_$m.fastq.gz.tmp" "${out}_$m.fastq.gz"
    done
}
