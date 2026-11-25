// hmm-annotate's Pfam domains of every MGnify protein -> the lookup bundle, a family group per task.
// See README.md.

process MEMBERSHIP {  // the release's cluster table is sorted by representative: here, by protein
    label 'process_host'

    input:
    path code

    output:
    path 'membership'

    script:
    """
    ${params.python} ${code} --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' \\
        membership --release '${params.release}' --shard-width ${params.shard_width}
    """

    stub:
    "mkdir membership"
}

process GROUPS {
    label 'process_host'
    publishDir params.outdir, mode: 'copy'

    input:
    path hmm, stageAs: 'hmm/*'
    path pfam, stageAs: 'pfam'
    path labels, stageAs: 'labels/*'
    val label_args
    path code

    output:
    path 'groups.tsv', emit: groups
    path 'batches.tsv', emit: batches

    script:
    """
    ${params.python} ${code} --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' groups \\
        --hmm ${hmm} --pfam pfam ${label_args} --groups ${params.groups} --batch-chunks ${params.batch_chunks}
    """

    stub:
    "printf 'model\\trows\\tgroup\\nPF00001\\t1\\t0\\n' > groups.tsv; printf 'batch\\tlo\\thi\\n0\\t0\\t1\\n' > batches.tsv"
}

process SPLIT {
    tag "${batch}"
    label 'process_split'

    input:
    tuple val(batch), val(lo), val(hi)
    path groups
    path membership
    path pfam, stageAs: 'pfam'
    path labels, stageAs: 'labels/*'
    val label_args
    path code

    output:
    path "batch${batch}"

    script:
    """
    ${params.python} ${code} --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' split \\
        --pfam pfam ${label_args} --groups ${groups} --membership ${membership} \\
        --lo ${lo} --hi ${hi} --out batch${batch}
    rm -rf duckdb_tmp
    """

    stub:
    "mkdir batch${batch}"
}

process BUILD {
    tag "${group}"
    label 'process_build'
    label 'lookup'

    input:
    val group
    path batches
    path hmm, stageAs: 'hmm/*'
    path code

    output:
    path "part${group}"

    script:
    """
    ${hmm_file(hmm)}
    ${params.lookup_python} ${code} build-part --rows batch*/group=${group}/*.parquet --hmm models.hmm --out part${group}
    """

    stub:
    "mkdir part${group}"
}

process ASSEMBLE {
    label 'process_assemble'
    label 'lookup'
    publishDir params.outdir, mode: params.publish_mode

    input:
    path parts
    path hmm, stageAs: 'hmm/*'
    path code

    output:
    path 'bundle'

    script:
    """
    ${hmm_file(hmm)}
    ${params.lookup_python} ${code} assemble --parts \$(ls -d part* | sort -V) --hmm models.hmm --out bundle \\
        --pfam '${params.pfam_version}' --mgnify '${params.mgnify_version}'
    rm models.hmm
    """

    stub:
    "mkdir bundle"
}

process CHECK {
    label 'process_build'
    label 'lookup'
    publishDir params.outdir, mode: 'copy'

    input:
    path bundle
    path sample
    path code

    output:
    path 'check*'

    script:
    def profile = sample ? "${params.lookup_python} ${code} profile --bundle ${bundle} --parquet ${sample} --out check --threads ${task.cpus}" : ''
    """
    ${params.lookup_python} ${code} verify ${bundle} > check.verify.txt
    ${profile}
    """

    stub:
    "touch check.verify.txt"
}

def hmm_file(hmm) {
    hmm.name.endsWith('.gz') ? "gunzip -c ${hmm} > models.hmm" : "cp -L ${hmm} models.hmm"
}

workflow {
    ['pfam', 'hmm'].each { p -> if (!params[p]) { error "--${p} is needed: see workflows/lookup-index/README.md" } }
    if (workflow.containerEngine && !params.lookup_sif) {
        error "--lookup_sif: the lookup image (workflows/pfam-anchor/lookup/Dockerfile)"
    }
    // ko=DIR,ec=DIR: hmm-annotate outputs whose models are the labels
    def specs = params.labels ? params.labels.toString().tokenize(',').collect { s -> s.tokenize('=') } : []
    if (specs.any { s -> s.size() != 2 }) {
        error "--labels ${params.labels}: scheme=directory, separated by commas"
    }
    def dirs = specs.collect { s -> file(s[1], checkIfExists: true) }
    if (dirs*.name.toSet().size() < dirs.size()) {
        error "--labels: the directories' names must differ"
    }
    def label_args = specs.indices.collect { i -> "--label ${specs[i][0]}=labels/${dirs[i].name}" }.join(' ')
    def hmm = file(params.hmm, checkIfExists: true)
    def pfam = file(params.pfam, checkIfExists: true)
    // inputs, so that containers see them and -resume reruns what an edit changes
    def code = file("${projectDir}/lookup_index.py", checkIfExists: true)
    def bundle_py = file("${projectDir}/../pfam-anchor/lookup/bundle.py", checkIfExists: true)
    def sample = params.check_parquet ? file(params.check_parquet, checkIfExists: true) : []

    MEMBERSHIP(file("${projectDir}/../mgnify-subset/mgnify_subset.py", checkIfExists: true))
    GROUPS(hmm, pfam, dirs, label_args, code)
    ch_batches = GROUPS.out.batches.splitCsv(sep: '\t', header: true).map { r -> [r.batch as int, r.lo as long, r.hi as long] }
    SPLIT(ch_batches, GROUPS.out.groups, MEMBERSHIP.out, pfam, dirs, label_args, code)
    ch_groups = GROUPS.out.groups.splitCsv(sep: '\t', header: true).map { r -> r.group as int }.unique()
    BUILD(ch_groups, SPLIT.out.collect(), hmm, bundle_py)
    ASSEMBLE(BUILD.out.collect(), hmm, bundle_py)
    CHECK(ASSEMBLE.out, sample, bundle_py)
}
