// aai-model: the parameters of the k-mer survival model aai inverts (plan, phase 7, steps 25
// and 33), fitted on MGnify protein pairs at nearest-member identity. Steps are in fit_model.py.

include { DB as REP_DB; DB as MEMBERS_DB; ALIGN as REP_ALIGN; ALIGN as NEAREST_ALIGN } from './diamond.nf'

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

process CANDIDATES {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'

    input:
    path hits

    output:
    path 'candidates.parquet'

    script:
    "${params.fit_model} candidates --hits ${hits} --per-bin ${params.per_bin_candidates} --min-cov ${params.min_cov} --seed ${params.seed}"

    stub:
    "touch candidates.parquet"
}

process MEMBER_DB {
    label 'process_medium'

    input:
    tuple path(candidates), path(members)

    output:
    path 'members.faa', emit: faa
    path 'member_clusters.parquet', emit: clusters
    path 'nearest_*.faa', emit: queries

    script:
    "${params.fit_model} member-db --candidates ${candidates} --members ${members} --chunks ${params.chunks}"

    stub:
    "touch members.faa member_clusters.parquet nearest_000.faa nearest_001.faa"
}

process PAIRS {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'

    input:
    path hits
    path candidates
    path member_clusters

    output:
    path 'pairs.parquet', emit: pairs
    path 'pairs_stats.json'

    script:
    """
    ${params.fit_model} pairs --hits ${hits} --candidates ${candidates} --member-clusters ${member_clusters} \\
        --per-bin ${params.per_bin} --min-cov ${params.min_cov} --seed ${params.seed}
    """

    stub:
    "touch pairs.parquet pairs_stats.json"
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
    // pass 1: P against the representatives finds candidate clusters at every identity
    REP_ALIGN(REP_DB(QUERIES.out.reps).combine(QUERIES.out.queries.flatten()), params.diamond_args)
    CANDIDATES(REP_ALIGN.out.collect())
    // pass 2: P against the candidates' members gives its nearest member of each
    MEMBER_DB(CANDIDATES.out.combine(Channel.value(members)))
    NEAREST_ALIGN(MEMBERS_DB(MEMBER_DB.out.faa).combine(MEMBER_DB.out.queries.flatten()), params.diamond_nearest_args)
    PAIRS(NEAREST_ALIGN.out.collect(), CANDIDATES.out, MEMBER_DB.out.clusters)
    FIT(SURVIVAL(PAIRS.out.pairs.combine(Channel.value(members))))
}
