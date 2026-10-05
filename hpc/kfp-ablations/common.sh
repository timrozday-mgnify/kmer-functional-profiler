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
