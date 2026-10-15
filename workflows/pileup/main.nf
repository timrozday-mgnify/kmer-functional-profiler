// Read recruitment and crude pileups downstream of a kfp profile (plan, phase 13; the
// nucleotide follow-up of phase 14, step 4). Per sample: kfp query -> per-sample reference
// (kfp-pileup reference) -> DIAMOND blastx (amino-acid index) or bwa-mem2 mem (nucleotide
// index), chosen by the index's molecule -> read posteriors, pileups and per-unit identity
// (kfp-pileup call). See README.md.

process INDEX {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'

    input:
    path targets
    val opts  // options as inputs, not params in the script: Nextflow does not hash params

    output:
    path 'index', emit: index

    script:
    "${params.kfp} index ${targets} index --targets ${opts}"

    stub:
    "mkdir index"
}

process QUERY {
    tag "${sample}"
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/${it}" }

    input:
    tuple val(sample), path(r1), path(r2)
    path index
    val opts  // e.g. '--stats stats.json' (query_stats) and query_args

    output:
    tuple val(sample), path('profile.tsv'), emit: profile
    path 'stats.json', optional: true  // with --query_stats: stage times, hits per pair

    script:
    "${params.kfp} query ${index} ${r1} ${r2} --out profile.tsv --aai ${opts}"

    stub:
    "touch profile.tsv"
}

process REFERENCE {
    tag "${sample}"
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/${it}" }

    input:
    tuple val(sample), path(profile)
    path index
    path store
    path withhold  // units left out of the reference (spike-benchmark), or NO_WITHHOLD
    val opts

    output:
    tuple val(sample), path('ref'), env('MOLECULE'), emit: ref

    script:
    def stores = store.collect { "--store ${it}" }.join(' ')
    def held = withhold.name == 'NO_WITHHOLD' ? '' : "--withhold ${withhold}"
    """
    ${params.kfp_pileup} reference ${profile} ${index} ref ${opts} ${stores} ${held}
    MOLECULE=\$(cat ref/molecule.txt)
    """

    stub:
    "mkdir ref && echo protein > ref/molecule.txt && MOLECULE=protein"
}

process DIAMOND {
    tag "${sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.2.8--he361c42_0'

    input:
    tuple val(sample), path(ref), path(r1), path(r2)
    val opts

    output:
    tuple val(sample), val('diamond'), path('r1.tsv'), path('r2.tsv'), emit: alignments

    script:
    // every hit within --top % of the read's best: the posteriors need the runners-up
    """
    diamond makedb --in ${ref}/ref.fa --db ref --threads ${task.cpus}
    for m in 1 2; do
        reads=\$([ \$m = 1 ] && echo ${r1} || echo ${r2})
        diamond blastx --db ref.dmnd --query \$reads --out r\$m.tsv --threads ${task.cpus} \\
            ${opts}
    done
    """

    stub:
    "touch r1.tsv r2.tsv"
}

process BWAMEM2 {
    tag "${sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he70b90d_8'

    input:
    tuple val(sample), path(ref), path(r1), path(r2)
    val opts

    output:
    tuple val(sample), val('bwa'), path('aln.sam'), emit: alignments

    script:
    // -a: every alignment, secondaries included, for the posteriors
    """
    bwa-mem2 index -p ref ${ref}/ref.fa
    bwa-mem2 mem -a -t ${task.cpus} ${opts} ref ${r1} ${r2} > aln.sam
    """

    stub:
    "touch aln.sam"
}

process CALL {
    tag "${sample}"
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/${it}" }

    input:
    tuple val(sample), val(aligner), path(alignments), path(ref), path(profile)
    val opts

    output:
    tuple val(sample), path('pileup'), emit: pileup

    script:
    """
    ${params.kfp_pileup} call ${ref} pileup ${alignments} --aligner ${aligner} \\
        --profile ${profile} ${opts}
    """

    stub:
    "mkdir pileup"
}

// Reusable by other pipelines (workflows/read-xref, spike-benchmark): reads (sample, r1, r2),
// the index, the store's members tables (a list, possibly empty) and a file of units to
// withhold from the reference (NO_WITHHOLD for none).
workflow PILEUP {
    take:
    reads
    index
    store
    withhold

    main:
    def profiles = QUERY(reads, index,
        "${params.query_stats ? '--stats stats.json ' : ''}${params.query_args}").profile
    def refs = REFERENCE(profiles, index, store, withhold,
        "--min-present ${params.min_present} --members ${params.members}").ref
    def to_align = refs.join(reads).branch { sample, ref, molecule, r1, r2 ->
        nucleotide: molecule == 'nucleotide'
            return tuple(sample, ref, r1, r2)
        protein: true
            return tuple(sample, ref, r1, r2)
    }
    def diamond = "--outfmt 6 ${params.diamond_fields} --top ${params.diamond_top} ${params.diamond_args}"
    def alignments = DIAMOND(to_align.protein, diamond).alignments
        .map { sample, aligner, r1, r2 -> tuple(sample, aligner, [r1, r2]) }
        .mix(BWAMEM2(to_align.nucleotide, params.bwa_args).alignments)
    def ref_dirs = refs.map { sample, ref, molecule -> tuple(sample, ref) }
    def pileups = CALL(alignments.join(ref_dirs).join(profiles),
        "--min-depth ${params.min_depth} --min-freq ${params.min_freq} --novel ${params.novel}").pileup

    emit:
    profile = profiles
    ref = ref_dirs
    pileup = pileups
}

workflow {
    if (!params.samples) {
        error 'Set --samples: a CSV with columns sample, r1, r2'
    }
    if (!params.index && !params.targets) {
        error 'Set --index (a kfp index) or --targets (a FASTA to build one from)'
    }
    def sheet = file(params.samples, checkIfExists: true)
    // reads paths may be relative to the sample sheet
    def reads = channel.fromPath(sheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample, file(sheet.parent.resolve(row.r1), checkIfExists: true),
                file(sheet.parent.resolve(row.r2), checkIfExists: true))
        }
    def index = params.index ? channel.value(file(params.index, checkIfExists: true))
        : INDEX(file(params.targets, checkIfExists: true), params.index_args).index
    def store = params.store ? channel.fromPath(params.store, checkIfExists: true).collect()
        : channel.value([])
    PILEUP(reads, index, store, file("${projectDir}/assets/NO_WITHHOLD"))
}
