// Tree index, the standard tier and the high-recall one (docs/tree-index-plan.md, Progress log:
// the standard index restored), built in parts. See README.md; on HPC,
// hpc/kfp-ablations/tree-index-build/run.sh. Run as: nextflow run workflows/tree-index/build.nf
//
// FIXTURE    synthetic members when --members is '' (the test profile)
// PARTITION  proteins -> parts of ~--part_proteins (k 11 storage components and equal
//            sequences kept whole) -> parts/part*.parquet, partition.json (partition.py)
// BUILD      per part and tier: build.py (a leaf per dense 8-mer set, dense 8-mer topology,
//            spread minimizers and the k 7 ladder diffed on the spanning tree; recall adds a
//            Dayhoff layer on every leaf)
// MERGE      per tier: the parts -> index-<tier>/ (merge.py), query.py's format
// TARGETS    --bench_targets proteins with an index position (recall.py targets)  (--bench false: skip)
// SIMULATE   one strain per target at 85-70% identity, BLOSUM62 substitutions (spike.py simulate)
// QUERY      per tier: the strains' reads, chance seeds counted on reversed reads (query.py --decoy)
// SCORE      per tier: candidate recall against the truth (recall.py score)
// SUMMARY    -> bench.tsv: per tier, postings, decoy bar and recall by identity (recall.py summary)

nextflow.enable.dsl = 2

params.members       = ''         // members parquet (absolute path or glob), '' = fixture
params.part_proteins = 2000000    // proteins per part
// name:build.py options, ';'-separated (split at each tier's first ':')
params.tiers         = 'standard:--mode standard;recall:--mode recall'
params.bench         = true
params.bench_targets = 2000
// sim4's settings (docs/tree-index-plan.md, Progress log: recall at 70-85% identity)
params.bench_args    = '--present 1 --identities 0.85 0.8 0.75 0.7 --substitution blosum62 --depth-mu 0.3 --depth-sigma 0.6'

def py() { "${params.venv}/bin/python" }

process FIXTURE {
    label 'process_low'

    output:
    path 'members.bucket*.parquet'

    script:
    "${py()} ${projectDir}/phase0.py fixture --families 40 --members 8 --buckets 1"
}

process PARTITION {
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy', pattern: 'partition.json'

    input:
    path members, stageAs: 'members/*'
    path code, stageAs: 'code/*'
    val part_proteins  // an input, so -resume reruns when it changes

    output:
    path 'part*.parquet', emit: parts
    path 'partition.json'

    script:
    "${py()} code/partition.py --members members/* --out . --part-proteins ${part_proteins}"
}

process BUILD {
    tag "${tier} ${part.baseName}"
    label 'process_medium'

    input:
    tuple val(tier), val(args), path(part)
    path code, stageAs: 'code/*'

    output:
    tuple val(tier), path("${tier}-${part.baseName}")

    script:
    "${py()} code/build.py --members ${part} --out ${tier}-${part.baseName} ${args}"
}

process MERGE {
    tag "${tier}"
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', saveAs: { "index-${tier}" }

    input:
    tuple val(tier), path(parts, stageAs: 'parts/*')
    path code, stageAs: 'code/*'

    output:
    path "index-${tier}"

    script:
    "${py()} code/merge.py --parts parts/* --out index-${tier}"
}

process TARGETS {
    label 'process_medium'

    input:
    path members, stageAs: 'members/*'
    path index, stageAs: 'index'
    path code, stageAs: 'code/*'

    output:
    path 'targets.fa'

    script:
    "${py()} code/recall.py targets --members members/* --index index --n ${params.bench_targets} --out targets.fa"
}

process SIMULATE {
    label 'process_low'

    input:
    path targets
    path code, stageAs: 'code/*'

    output:
    path 'sim'

    script:
    "${py()} code/spike.py simulate --targets ${targets} ${params.bench_args} --out sim"
}

process QUERY {
    tag "${tier}"
    label 'process_high_memory'

    input:
    tuple val(tier), path(index, stageAs: 'index'), path(sim)
    path code, stageAs: 'code/*'

    output:
    tuple val(tier), path("q-${tier}")

    script:
    "${py()} code/query.py --index index --reads ${sim}/reads_1.fq.gz ${sim}/reads_2.fq.gz --decoy --out q-${tier}"
}

process SCORE {
    tag "${tier}"
    label 'process_medium'

    input:
    tuple val(tier), path(index, stageAs: 'index'), path(q), path(sim)
    path code, stageAs: 'code/*'

    output:
    tuple val(tier), path("s-${tier}")

    script:
    """
    ${py()} code/recall.py score --index index --query ${q} --truth ${sim}/truth.tsv \
        --out s-${tier} --caps 256 --deltas 0.1
    """
}

process SUMMARY {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path runs, stageAs: 'runs/*'
    path code, stageAs: 'code/*'

    output:
    path 'bench.tsv'
    path 'runs/[qs]-*', includeInputs: true

    script:
    "${py()} code/recall.py summary --root runs --base standard --out bench.tsv"
}

workflow {
    code = ['partition.py', 'build.py', 'merge.py', 'survey.py', 'phase0.py', 'query.py',
            'alphabets.py', 'recall.py', 'report.py', '../spike-benchmark/spike.py']
        .collect { file("${projectDir}/${it}") }
    members = params.members ? Channel.fromPath(params.members, checkIfExists: true).collect()
                             : FIXTURE()
    PARTITION(members, code, params.part_proteins)
    tiers = Channel.fromList(params.tiers.tokenize(';')*.trim().collect { it.split(':', 2) as List })
    BUILD(tiers.combine(PARTITION.out.parts.flatten()), code)
    MERGE(BUILD.out.groupTuple(), code)
    if (params.bench) {
        index = MERGE.out.map { [it.name - 'index-', it] }
        TARGETS(members, index.first().map { it[1] }, code)
        sim = SIMULATE(TARGETS.out, code)
        QUERY(index.combine(sim), code)
        SCORE(index.join(QUERY.out).combine(sim), code)
        // index-<tier>, q-<tier>, s-<tier> side by side
        SUMMARY(MERGE.out.mix(QUERY.out.map { it[1] }, SCORE.out.map { it[1] }).collect(), code)
    }
}
