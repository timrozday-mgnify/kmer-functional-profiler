// aai-model: the parameters of the k-mer survival model aai inverts (plan, phase 7, step 25),
// fitted on MGnify protein pairs. Steps are in fit_model.py.

process QUERIES {
    label 'process_medium'

    input:
    path members

    output:
    path 'reps.faa', emit: reps
    path 'queries_*.faa', emit: queries

    script:
    "${params.fit_model} queries --members ${members} --n ${params.n_queries} --chunks ${params.chunks} --seed ${params.seed}"

    stub:
    "touch reps.faa queries_000.faa queries_001.faa"
}

process DB {
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.2.8--he361c42_0'

    input:
    path faa

    output:
    path 'reps.dmnd'

    script:
    "diamond makedb --in ${faa} --db reps --threads ${task.cpus}"

    stub:
    "touch reps.dmnd"
}

process ALIGN {
    tag "${faa.baseName}"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.2.8--he361c42_0'

    input:
    tuple path(db), path(faa)

    output:
    path "${faa.baseName}.hits.tsv"

    script:
    """
    diamond blastp --db ${db} --query ${faa} --out ${faa.baseName}.hits.tsv --threads ${task.cpus} \\
        --outfmt 6 qseqid sseqid pident length qlen slen bitscore ${params.diamond_args}
    """

    stub:
    "touch ${faa.baseName}.hits.tsv"
}

process PAIRS {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'

    input:
    path hits

    output:
    path 'pairs.parquet'

    script:
    "${params.fit_model} pairs --hits ${hits} --per-bin ${params.per_bin} --min-cov ${params.min_cov} --seed ${params.seed}"

    stub:
    "touch pairs.parquet"
}

process SURVIVAL {
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(pairs), path(members)

    output:
    path 'survival.parquet'

    script:
    "${params.fit_model} survival --pairs ${pairs} --members ${members} --k ${params.k} --alphabet ${params.alphabet}"

    stub:
    "touch survival.parquet"
}

process FIT {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'

    input:
    path survival

    output:
    path 'aai_model.json'
    path 'model_strata.tsv'

    script:
    "${params.fit_model} fit --survival ${survival} --k ${params.k} --alphabet ${params.alphabet}"

    stub:
    "touch aai_model.json model_strata.tsv"
}

workflow {
    def members = file(params.members, checkIfExists: true)
    QUERIES(members)
    ALIGN(DB(QUERIES.out.reps).combine(QUERIES.out.queries.flatten()))
    PAIRS(ALIGN.out.collect())
    FIT(SURVIVAL(PAIRS.out.combine(Channel.value(members))))
}
