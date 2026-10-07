// DIAMOND database and blastp of the MGnify90-level truth (main.nf): once against the
// representatives, once against the hit clusters' members (nearest-member identity).

process MGNIFY_DB {
    tag "${name}"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.2.8--he361c42_0'

    input:
    tuple val(name), path(faa)

    output:
    tuple val(name), path("${faa.baseName}_${name}.dmnd"), emit: db

    script:
    "diamond makedb --in ${faa} --db ${faa.baseName}_${name} --threads ${task.cpus}"

    stub:
    "touch ${faa.baseName}_${name}.dmnd"
}

process MGNIFY_ANNOTATE {
    tag "${name} ${faa.baseName}"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.2.8--he361c42_0'

    input:
    tuple val(name), path(db), path(faa)
    val args

    output:
    tuple val(name), path("${faa.baseName}.hits.tsv"), emit: hits

    script:
    """
    diamond blastp --db ${db} --query ${faa} --out ${faa.baseName}.hits.tsv --threads ${task.cpus} \\
        --outfmt 6 qseqid sseqid pident length qlen slen bitscore ${args}
    """

    stub:
    "touch ${faa.baseName}.hits.tsv"
}
