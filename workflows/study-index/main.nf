// Study index (plan: Additional references, step 4): a study's contigs or MAGs -> pyrodigal
// (meta) -> MMseqs2 linclust at 90% -> members, Pfam by hmmsearch --cut_ga -> `index --like
// BASE`, to be queried jointly with the base (`query BASE READS --extra-index INDEX`).
// STUDY_INDEX is also included by the fmh benchmark's study ladder. See README.md.


process STUDY_GENES {
    tag "${name}"
    label 'process_medium'
    container 'quay.io/biocontainers/pyrodigal:3.7.1--py312h247cb63_1'

    input:
    tuple val(name), path(fastas, stageAs: 'in/*')

    output:
    tuple val(name), path("${name}.faa"), emit: faa

    script:
    // ids are prefixed with the file's name, so contig names may repeat across MAGs
    """
    for f in in/*; do
        stem=\$(basename "\$f"); stem=\${stem%%.*}
        pyrodigal -p meta -j ${task.cpus} -i "\$f" -a genes.faa -o /dev/null
        sed "s/^>/>\${stem}_/" genes.faa >> ${name}.faa
    done
    rm -f genes.faa
    """

    stub:
    "touch ${name}.faa"
}

process STUDY_PROTEINS {
    tag "${name}"
    label 'process_single'

    input:
    tuple val(name), path(fastas, stageAs: 'in/*')

    output:
    tuple val(name), path("${name}.faa"), emit: faa

    script:
    "for f in in/*; do gzip -cdf \"\$f\"; done > ${name}.faa"

    stub:
    "touch ${name}.faa"
}

process STUDY_CLUSTER {
    tag "${name}"
    label 'process_medium'
    container 'quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0'

    input:
    tuple val(name), path(faa)
    val args

    output:
    tuple val(name), path("${name}.clusters.tsv"), emit: clusters

    script:
    """
    mmseqs createdb ${faa} db
    mmseqs linclust db clu tmp ${args} --threads ${task.cpus}
    mmseqs createtsv db db clu ${name}.clusters.tsv
    rm -rf tmp db* clu*
    """

    stub:
    "touch ${name}.clusters.tsv"
}

process STUDY_MEMBERS {
    tag "${name}"
    label 'process_single'

    input:
    tuple val(name), path(faa), path(clusters)

    output:
    tuple val(name), path("${name}.members.parquet"), emit: members

    script:
    "${params.python} ${moduleDir}/study.py members --faa ${faa} --clusters ${clusters} --out ${name}.members.parquet"

    stub:
    "touch ${name}.members.parquet"
}

process STUDY_PFAM_DB {
    label 'process_single'
    storeDir "${params.db_dir}/pfam"

    output:
    path 'Pfam-A.hmm', emit: hmm

    script:
    "curl -fsSL ${params.pfam_url}/Pfam-A.hmm.gz | gunzip > Pfam-A.hmm"

    stub:
    "touch Pfam-A.hmm"
}

process STUDY_PFAM_SEARCH {
    tag "${name}:${chunk}"
    label 'process_medium'
    container 'quay.io/biocontainers/hmmer:3.4--hdbdd923_2'

    input:
    tuple val(name), val(chunk), val(chunks), path(faa)
    path hmm

    output:
    tuple val(name), path("${name}.${chunk}.domtbl"), emit: domtbl

    script:
    // every chunks-th protein from the chunk-th, so jobs get equal shares
    """
    awk -v n=${chunks} -v c=${chunk} '/^>/ { keep = (i++ % n == c) } keep' ${faa} > part.faa
    hmmsearch --cpu ${task.cpus} --cut_ga --domtblout ${name}.${chunk}.domtbl -o /dev/null ${hmm} part.faa
    rm part.faa
    """

    stub:
    "touch ${name}.${chunk}.domtbl"
}

process STUDY_PFAM {
    tag "${name}"
    label 'process_single'

    input:
    tuple val(name), path(domtbl, stageAs: 'tbl/*')

    output:
    tuple val(name), path("${name}.pfam.parquet"), emit: pfam

    script:
    "${params.python} ${moduleDir}/study.py pfam --domtbl tbl/* --out ${name}.pfam.parquet"

    stub:
    "touch ${name}.pfam.parquet"
}

process STUDY_INDEX_BUILD {
    tag "${name}"
    label 'process_index'
    publishDir "${params.outdir}/study", mode: 'copy'

    input:
    tuple val(name), path(members), path(pfam), path(base)
    val args

    output:
    tuple val(name), path("${name}"), emit: index
    tuple val(name), path("${name}.stats.json"), emit: stats

    script:
    "${params.kfp} index ${members} ${name} --like ${base} --pfam ${pfam} ${args} > ${name}.stats.json"

    stub:
    "mkdir ${name} && touch ${name}/meta.json ${name}.stats.json"
}

workflow STUDY_INDEX {
    take:
    contigs       // [name, [nucleotide FASTA]]: contigs or MAGs, plain or gzip
    proteins      // [name, [protein FASTA]]: gene calls already made
    base          // index directory whose build parameters the study index copies
    hmm           // Pfam-A.hmm
    cluster_args  // mmseqs linclust options
    pfam_chunks   // hmmsearch jobs per study (a value channel)
    index_args    // more `index` options, e.g. '--role decoy'

    main:
    faa = STUDY_GENES(contigs).faa.mix(STUDY_PROTEINS(proteins).faa)
    clusters = STUDY_CLUSTER(faa, cluster_args).clusters
    members = STUDY_MEMBERS(faa.join(clusters)).members
    chunks = faa.combine(pfam_chunks).flatMap { name, f, n -> (0..<n).collect { c -> [name, c, n, f] } }
    domtbl = STUDY_PFAM_SEARCH(chunks, hmm).domtbl
    pfam = STUDY_PFAM(domtbl.groupTuple()).pfam
    built = STUDY_INDEX_BUILD(
        members.join(pfam).combine(base),
        index_args,
    )

    emit:
    index = built.index
    stats = built.stats
}

workflow {
    if (!params.base) error 'set --base, the index the study index is queried with'
    if (!params.contigs && !params.proteins) error 'set --contigs or --proteins'
    def contigs = params.contigs ? Channel.of([params.name, files(params.contigs, checkIfExists: true)]) : Channel.empty()
    def proteins = params.proteins ? Channel.of([params.name, files(params.proteins, checkIfExists: true)]) : Channel.empty()
    def hmm = params.pfam_hmm ? Channel.value(file(params.pfam_hmm, checkIfExists: true)) : STUDY_PFAM_DB().hmm
    STUDY_INDEX(
        contigs,
        proteins,
        Channel.value(file(params.base, checkIfExists: true)),
        hmm,
        params.cluster_args,
        Channel.value(params.pfam_chunks as int),
        params.role ? "--role ${params.role}" : '',
    )
}
