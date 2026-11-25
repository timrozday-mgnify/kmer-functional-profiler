// nail's seeds (an MMseqs2 profile search) as candidate (model, protein) pairs, for annotate.py --seeds.
// A module so that it can be included twice: representatives, and their members.

process NAIL_SEED {
    tag "${chunk}"
    label 'process_nail'
    publishDir "${params.outdir}/logs", mode: params.publish_mode, pattern: '*.log'

    input:
    tuple val(chunk), path(faa)
    path hmm

    output:
    tuple val(chunk), path("chunk${chunk}.seeds.*.tsv.gz"), emit: seeds
    path "nail${chunk}.log", emit: log

    // nail's temporary files and memory grow with the proteins, and it spends a minute or more on the models
    // whatever they are: --nail_proteins at a time. `time` counts the mmseqs processes nail runs.
    script:
    """
    set -o pipefail
    trap 'rm -rf tmp part.*; cat nail${chunk}.log >&2' EXIT
    TIMEFORMAT='cpu_s %U %S'
    { grep -m1 'model name' /proc/cpuinfo || uname -m; } > nail${chunk}.log 2> /dev/null
    split -a 3 -l ${2 * (params.nail_proteins as long)} ${faa} part.
    for part in part.*; do
        { time ${params.nail} search --only-seed -S ${params.nail_seed_pvalue} -t ${task.cpus} --tmp-dir tmp \\
            --seeds-out seeds.tsv models.hmm \$part; } >> nail${chunk}.log 2>&1
        { printf 'model\\tprotein_id\\n'; cut -f 1,2 seeds.tsv; } | gzip -1 > chunk${chunk}.seeds.\${part#part.}.tsv.gz
        rm -rf tmp seeds.tsv
    done
    """

    stub:
    "touch nail${chunk}.log; echo | gzip > chunk${chunk}.seeds.aaa.tsv.gz"
}
