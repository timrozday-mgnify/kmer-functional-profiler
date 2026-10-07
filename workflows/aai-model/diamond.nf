// DIAMOND database and blastp, used for both of aai-model's passes (main.nf).

process DB {
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.2.8--he361c42_0'

    input:
    path faa

    output:
    path "${faa.baseName}.dmnd"

    script:
    "diamond makedb --in ${faa} --db ${faa.baseName} --threads ${task.cpus}"

    stub:
    "touch ${faa.baseName}.dmnd"
}

process ALIGN {
    tag "${faa.baseName}"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.2.8--he361c42_0'

    input:
    tuple path(db), path(faa)
    val args

    output:
    path "${faa.baseName}.hits.tsv"

    script:
    """
    diamond blastp --db ${db} --query ${faa} --out ${faa.baseName}.hits.tsv --threads ${task.cpus} \\
        --outfmt 6 qseqid sseqid pident length qlen slen bitscore ${args}
    """

    stub:
    "touch ${faa.baseName}.hits.tsv"
}
