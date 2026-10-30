// Tree index, phase 0: measure the k-mer tree on a sample of MGnify, no index
// (docs/tree-index-plan.md, Implementation plan, phase 0). See README.md; on HPC,
// hpc/kfp-ablations/tree-index-phase0/run.sh.
//
// SAMPLE   per bucket: clusters with cluster_rep % --sample == 0, their stored k-mers
// CORE     -> sample.parquet, core_hashes.npy
// STREAM   per bucket of the whole release: set digests, signatures, k-mer slices by hash
//          range; representatives outside the sample sharing stored k-mers with it
// REDUCE   per hash range: group sizes                                (--stream false: skip)
// COMBINE  -> stream.tsv, stream.json
// EXTRA    -> extra.parquet (the representatives)
// MEASURE  the sample (core), and the sample plus representatives (all) -> measure-*/
// TIER1    today's index on the sample's proteins -> tier1/meta.json
// REPORT   -> gate.tsv
// TIER1_MATCHED  today's index at the survey's other rates          (--survey false: skip)
// SURVEY   the floor survey on core and all -> survey-*/ (survey.py)
// INDEX_SIZES    file sizes of the release's standard tier (--release_index)

nextflow.enable.dsl = 2

def phase0() {
    "${params.venv}/bin/python ${projectDir}/phase0.py"
}

process FIXTURE {
    label 'process_low'

    output:
    path 'members.bucket*.parquet'

    script:
    "${phase0()} fixture ${params.tree_args}"
}

process SAMPLE {
    tag "${members.baseName}"
    label 'process_low'

    input:
    path members
    path code

    output:
    path "${members.baseName}.members.parquet", emit: members
    path "${members.baseName}.hashes.npy", emit: hashes

    script:
    """
    ${phase0()} sample --members ${members} --sample ${params.sample} \\
        --prefix ${members.baseName} ${params.tree_args}
    """
}

process CORE {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path members
    path hashes

    output:
    path 'sample.parquet', emit: members
    path 'core_hashes.npy', emit: hashes

    script:
    """
    ${phase0()} concat --members ${members} --hashes ${hashes} --out sample.parquet \\
        --out-hashes core_hashes.npy
    """
}

process STREAM {
    tag "${members.baseName}"
    label 'process_low'

    input:
    path members
    path core
    path code

    output:
    path "${members.baseName}.r*.parquet", emit: parts
    path "${members.baseName}.extra.parquet", emit: extra
    path "${members.baseName}.totals.json", emit: totals

    script:
    """
    ${phase0()} stream --members ${members} --prefix ${members.baseName} --core-hashes ${core} \\
        --sample ${params.sample} --ranges ${params.ranges} ${params.tree_args} ${params.stream_args}
    """
}

process REDUCE {
    tag "range ${range}"
    label 'process_medium'

    input:
    tuple val(range), path(parts)
    val passes  // an input, so -resume reruns when it changes (params are not hashed)

    output:
    path "r${range}.hist.parquet", emit: hist
    path "r${range}.top.parquet", emit: top

    script:  // reduce.py, not phase0.py: editing it keeps SAMPLE and STREAM cached
    """
    ${params.venv}/bin/python ${projectDir}/reduce.py --parts ${parts} --prefix r${range} \
        --passes ${passes}
    """
}

process COMBINE {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path hist
    path top
    path totals

    output:
    path 'stream.tsv'
    path 'stream.json'

    script:
    "${phase0()} combine --hist ${hist} --tops ${top} --totals ${totals}"
}

process EXTRA {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path extra

    output:
    path 'extra.parquet'

    script:
    "${phase0()} concat --members ${extra} --out extra.parquet --max-rows ${params.max_extra}"
}

process MEASURE {
    tag "${name}"
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(name), path(members), path(extra)
    path code

    output:
    path "measure-${name}", emit: out

    script:
    def more = extra ? "--extra ${extra}" : ''
    """
    ${phase0()} measure --members ${members} ${more} --out measure-${name} \\
        ${params.tree_args} ${params.measure_args}
    """
}

process TIER1 {
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy', pattern: 'tier1/meta.json'

    input:
    path members
    val index_args  // an input, so -resume rebuilds when it changes

    output:
    path 'tier1/meta.json', emit: meta

    script:
    "${params.kfp} index ${members} tier1 ${index_args}"
}

process TIER1_MATCHED {
    tag "${index_args}"
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy', saveAs: { "tier1-${index_args.replaceAll(/[^0-9a-z.]+/, '_')}.json" }

    input:
    path members
    each index_args

    output:
    path 'tier1/meta.json', emit: meta

    script:
    "${params.kfp} index ${members} tier1 ${index_args}"
}

process SURVEY {
    tag "${name}"
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(name), path(members), path(extra)
    path metas, stageAs: 'tier1-*.json'
    path stream
    path code, stageAs: 'code/*'
    val survey_args  // an input, so -resume reruns when it changes

    output:
    path "survey-${name}"

    script:
    def more = extra ? "--extra ${extra}" : ''
    def st = stream ? "--stream ${stream}" : ''
    """
    ${params.venv}/bin/python code/survey.py --members ${members} ${more} --out survey-${name} \
        --tier1 ${metas} ${st} ${survey_args}
    """
}

process INDEX_SIZES {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    val index

    output:
    path 'index_sizes.tsv'

    script:
    """
    printf 'file\tbytes\n' > index_sizes.tsv
    find -L ${index} -maxdepth 1 -type f -printf '%f\t%s\n' | sort >> index_sizes.tsv
    """
}

process REPORT {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path measure
    path meta
    path code

    output:
    path 'gate.tsv'

    script:
    "${phase0()} report --measure ${measure}/measure.json --tier1 ${meta} --out gate.tsv"
}

workflow {
    code = file("${projectDir}/phase0.py")
    buckets = params.store ? Channel.fromPath(params.store, checkIfExists: true) : FIXTURE().flatten()
    SAMPLE(buckets, code)
    CORE(SAMPLE.out.members.collect(), SAMPLE.out.hashes.collect())
    extra = Channel.value([])
    if (params.stream) {
        STREAM(buckets, CORE.out.hashes, code)
        by_range = STREAM.out.parts.flatten()
            .map { f -> [(f.name =~ /\.r(\d+)\.parquet$/)[0][1] as int, f] }
            .groupTuple()
        REDUCE(by_range, params.reduce_passes)
        COMBINE(REDUCE.out.hist.collect(), REDUCE.out.top.collect(), STREAM.out.totals.collect())
        extra = EXTRA(STREAM.out.extra.collect())
    }
    runs = CORE.out.members.map { m -> ['core', m, []] }
    if (params.stream) {
        runs = runs.mix(CORE.out.members.combine(extra).map { m, e -> ['all', m, e] })
    }
    MEASURE(runs, code)
    TIER1(CORE.out.members, params.index_args)
    core = MEASURE.out.out.filter { it.name == 'measure-core' }
    REPORT(core, TIER1.out.meta, code)
    if (params.survey) {
        TIER1_MATCHED(CORE.out.members, params.tier1_variants.tokenize(';')*.trim())
        metas = TIER1.out.meta.mix(TIER1_MATCHED.out.meta).collect()
        stream = params.stream ? COMBINE.out[1] : Channel.value([])
        survey_code = [file("${projectDir}/survey.py"), file("${projectDir}/phase0.py")]
        SURVEY(runs, metas, stream, survey_code, params.survey_args)
    }
    if (params.release_index) {
        INDEX_SIZES(params.release_index)
    }
}
