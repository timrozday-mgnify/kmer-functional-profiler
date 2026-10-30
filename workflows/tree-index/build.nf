// Tree index, the three tiers (docs/tree-index-plan.md, Progress log: tiers chosen from the
// per-protein survey), built in parts. See README.md; on HPC,
// hpc/kfp-ablations/tree-index-build/run.sh. Run as: nextflow run workflows/tree-index/build.nf
//
// FIXTURE    synthetic members when --members is '' (the test profile)
// PARTITION  proteins -> parts of ~--part_proteins (k 11 storage components and equal
//            sequences kept whole) -> parts/part*.parquet, partition.json (partition.py)
// BUILD      per part and tier: build.py (a leaf per dense 8-mer set, dense 8-mer topology)
// MERGE      per tier: the parts -> index-<tier>/ (merge.py), query.py's format

nextflow.enable.dsl = 2

params.members       = ''         // members parquet (absolute path or glob), '' = fixture
params.part_proteins = 2000000    // proteins per part
// name:build.py options, ';'-separated
params.tiers         = 'quick:--t-base 0.01 --n-min 8 --tau 0 --n-1 1;' +
                       'medium:--t-base 0.02 --n-min 16 --tau 0.3 --n-1 2;' +
                       'large:--t-base 0.05 --n-min 32 --tau 0.3 --n-1 2'

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

workflow {
    code = ['partition.py', 'build.py', 'merge.py', 'survey.py', 'phase0.py']
        .collect { file("${projectDir}/${it}") }
    members = params.members ? Channel.fromPath(params.members, checkIfExists: true).collect()
                             : FIXTURE()
    PARTITION(members, code, params.part_proteins)
    tiers = Channel.fromList(params.tiers.tokenize(';')*.trim().collect { it.tokenize(':') })
    BUILD(tiers.combine(PARTITION.out.parts.flatten()), code)
    MERGE(BUILD.out.groupTuple(), code)
}
