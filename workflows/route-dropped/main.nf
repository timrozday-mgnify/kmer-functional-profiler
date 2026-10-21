// Routing's dropped reads against every profile unit (plan, phase 13, routing): DROPPED
// (the reads, single-end, routing does not send to DIAMOND, subsampled), FULLREF (every profile unit,
// 1 and 3 members), DIAMOND (as workflows/pileup's: --sensitive -F 15 --top 10) of the
// dropped reads against each full reference. dropped.py summarise reads the outputs. See
// README.md; on HPC, hpc/kfp-ablations/route-dropped/run.sh.

// The venv's python; an empty params.python (a stale config setting it to '') falls back to
// it rather than running dropped.py as a command
def python() {
    params.python ?: "${params.venv}/bin/python"
}

process DROPPED {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'

    input:
    path ref
    path meta, stageAs: 'index/meta.json'
    tuple path(r1), path(r2)
    path code

    output:
    path 'dropped', emit: dir

    script:
    "${python()} ${code} dropped --ref ${ref} --index index --reads ${r1} ${r2} --reads-n ${params.reads} --out dropped"
}

process FULLREF {
    tag "m${members}"
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', pattern: '*.tsv', saveAs: { "m${members}/${it}" }

    input:
    path profile
    path ref
    val index  // the index directory, read in place (its units table and store)
    val members
    path code

    output:
    tuple val("m${members}"), path('full.fa'), path('full.tsv'), emit: ref

    script:
    def store = params.store ? "--store ${params.store}" : ''
    def held = params.withhold ? "--withhold ${params.withhold}" : ''
    "${python()} ${code} fullref --profile ${profile} --ref ${ref} --index ${index} --members ${members} ${store} ${held}"
}

process DIAMOND {
    tag "${arm}"
    label 'process_high'
    container 'quay.io/biocontainers/diamond:2.2.8--he361c42_0'
    publishDir params.outdir, mode: 'copy', pattern: '*.tsv', saveAs: { "${arm}/${it}" }

    input:
    tuple val(arm), path(fasta), path(table), path(dropped)

    output:
    tuple val(arm), path('r.tsv'), path(table)

    script:
    // dropped reads are single-end (each mate routed on its own)
    """
    diamond makedb --in ${fasta} --db full --threads ${task.cpus}
    diamond blastx --db full.dmnd --query ${dropped}/dropped.fa.gz --out r.tsv \\
        --threads ${task.cpus} --outfmt 6 qseqid sseqid pident length bitscore \\
        --top 10 -F 15 --sensitive
    """
}

workflow {
    def code = file("${projectDir}/dropped.py", checkIfExists: true)
    def reads = tuple(file(params.r1, checkIfExists: true), file(params.r2, checkIfExists: true))
    def ref = file(params.ref, checkIfExists: true)
    def dropped = DROPPED(ref, file("${params.index}/meta.json", checkIfExists: true), reads, code).dir
    def refs = FULLREF(file(params.profile, checkIfExists: true), ref, params.index,
        channel.of(1, 3), code).ref
    DIAMOND(refs.combine(dropped))
}
