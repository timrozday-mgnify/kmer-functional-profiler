// MGnify proteins -> hmmsearch against Pfam-A -> every domain hit, marked above or below GA.
// See README.md.

process PFAM_HMM {
    label 'process_small'
    publishDir params.outdir, mode: 'copy', pattern: 'ga.tsv'

    input:
    path hmm

    output:
    tuple path('Pfam-A.hmm'), path('ga.tsv')

    script:
    """
    if [[ ${hmm} == *.gz ]]; then gunzip -c ${hmm} > tmp.hmm; else cp -L ${hmm} tmp.hmm; fi
    mv tmp.hmm Pfam-A.hmm
    ${params.python} ${projectDir}/reannotate.py ga --hmm Pfam-A.hmm
    """

    stub:
    "touch Pfam-A.hmm ga.tsv"
}

process PLAN {
    label 'process_plan'
    publishDir params.outdir, mode: 'copy'

    output:
    path 'ranges.tsv'

    script:
    """
    ${params.python} ${projectDir}/reannotate.py --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' plan --release '${params.release}' \\
        --targets ${params.targets} --sample ${params.sample} \\
        --chunk-proteins ${params.chunk_proteins} ${params.plan_args}
    """

    stub:
    "printf 'chunk\\tlo\\thi\\tproteins\\n0\\t0\\t100\\t1\\n' > ranges.tsv"
}

process FASTA {
    tag "${arm} ${chunk}"
    label 'process_small'

    input:
    tuple val(arm), val(chunk), val(lo), val(hi), val(proteins)

    output:
    tuple val(arm), val(chunk), path("${arm}.chunk${chunk}.faa"), path("${arm}.chunk${chunk}.faa.tsv")

    script:
    """
    ${params.python} ${projectDir}/reannotate.py --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' fasta --release '${params.release}' \\
        --targets ${params.targets} --sample ${params.sample} ${arm == 'reversed' ? '--reverse' : ''} \\
        --lo ${lo} --hi ${hi} --proteins ${proteins} --out ${arm}.chunk${chunk}.faa
    """

    stub:
    "touch ${arm}.chunk${chunk}.faa ${arm}.chunk${chunk}.faa.tsv"
}

process HMMSEARCH {
    tag "${arm} ${chunk}"
    label 'process_hmmsearch'
    container 'quay.io/biocontainers/hmmer:3.4--hdbdd923_2'

    input:
    tuple val(arm), val(chunk), path(faa), path(counts)
    tuple path(hmm), path(ga)

    output:
    tuple val(arm), val(chunk), path("${arm}.chunk${chunk}.domtbl.gz"), path("${arm}.chunk${chunk}.txt.gz"), path(counts)

    script:
    """
    set -o pipefail
    hmmsearch --notextw --cpu ${task.cpus} ${params.hmmsearch_args} --domtblout tmp.domtbl ${hmm} ${faa} | gzip -1 > tmp.txt.gz
    gzip -1 tmp.domtbl
    mv tmp.domtbl.gz ${arm}.chunk${chunk}.domtbl.gz
    mv tmp.txt.gz ${arm}.chunk${chunk}.txt.gz
    """

    stub:
    "touch ${arm}.chunk${chunk}.domtbl.gz ${arm}.chunk${chunk}.txt.gz"
}

process PARSE {
    tag "${arm} ${chunk}"
    label 'process_small'
    publishDir "${params.outdir}/hits", mode: params.publish_mode, pattern: '*.parquet'

    input:
    tuple val(arm), val(chunk), path(domtbl), path(alignments), path(counts)
    tuple path(hmm), path(ga)

    output:
    path "${arm}.chunk${chunk}.parquet", emit: hits
    path "${arm}.chunk${chunk}.stats.tsv", emit: stats

    script:
    """
    ${params.python} ${projectDir}/reannotate.py parse --arm ${arm} --domtbl ${domtbl} --alignments ${alignments} \\
        ${flag(params.pp) ? '--pp' : ''} --ga ${ga} --counts ${counts} \\
        --out ${arm}.chunk${chunk}.parquet --stats ${arm}.chunk${chunk}.stats.tsv
    """

    stub:
    "touch ${arm}.chunk${chunk}.parquet ${arm}.chunk${chunk}.stats.tsv"
}

process SUMMARY {
    label 'process_small'
    publishDir params.outdir, mode: 'copy'

    input:
    path stats

    output:
    path 'summary.tsv'

    script:
    "${params.python} ${projectDir}/reannotate.py summary ${stats}"

    stub:
    "touch summary.tsv"
}

process COMPARE {
    label 'process_small'
    publishDir params.outdir, mode: 'copy'

    input:
    path hits

    output:
    path 'compare.tsv'

    script:
    "${params.python} ${projectDir}/reannotate.py compare ${hits}"

    stub:
    "touch compare.tsv"
}

// A boolean param: Nextflow 26 keeps `--reversed false` from the command line as the string
// "false", which Groovy reads as true; a config's `false` stays a Boolean.
def flag(value) {
    value instanceof String ? value.toBoolean() : value as boolean
}

workflow {
    if (!(params.targets in ['reps', 'all'])) {
        error "--targets ${params.targets}: reps or all"
    }
    PFAM_HMM(channel.value(file(params.hmm, checkIfExists: true)))
    ch_hmm = PFAM_HMM.out.first()
    ch_chunks = PLAN().splitCsv(sep: '\t', header: true)
        .map { r -> [r.chunk as int, r.lo as long, r.hi as long, r.proteins as long] }
    def arms = flag(params.reversed) ? ['real', 'reversed'] : ['real']
    FASTA(channel.fromList(arms).combine(ch_chunks))
    HMMSEARCH(FASTA.out, ch_hmm)
    PARSE(HMMSEARCH.out, ch_hmm)
    SUMMARY(PARSE.out.stats.collect())
    if (flag(params.reversed)) {
        COMPARE(PARSE.out.hits.collect())
    }
}
