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

    output:
    path 'index', emit: index

    script:
    "${params.kfp} index ${targets} index --targets ${params.index_args}"

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

    output:
    tuple val(sample), path('profile.tsv'), emit: profile

    script:
    "${params.kfp} query ${index} ${r1} ${r2} --out profile.tsv --aai ${params.query_args}"

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

    output:
    tuple val(sample), path('ref'), env('MOLECULE'), emit: ref

    script:
    def stores = store.collect { "--store ${it}" }.join(' ')
    """
    ${params.kfp_pileup} reference ${profile} ${index} ref --min-present ${params.min_present} \\
        --members ${params.members} ${stores}
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

    output:
    tuple val(sample), val('diamond'), path('r1.tsv'), path('r2.tsv'), emit: alignments

    script:
    // every hit within --top % of the read's best: the posteriors need the runners-up
    """
    diamond makedb --in ${ref}/ref.fa --db ref --threads ${task.cpus}
    for m in 1 2; do
        reads=\$([ \$m = 1 ] && echo ${r1} || echo ${r2})
        diamond blastx --db ref.dmnd --query \$reads --out r\$m.tsv --threads ${task.cpus} \\
            --outfmt 6 ${params.diamond_fields} --top ${params.diamond_top} ${params.diamond_args}
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

    output:
    tuple val(sample), val('bwa'), path('aln.sam'), emit: alignments

    script:
    // -a: every alignment, secondaries included, for the posteriors
    """
    bwa-mem2 index -p ref ${ref}/ref.fa
    bwa-mem2 mem -a -t ${task.cpus} ${params.bwa_args} ref ${r1} ${r2} > aln.sam
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

    output:
    tuple val(sample), path('pileup'), emit: pileup

    script:
    """
    ${params.kfp_pileup} call ${ref} pileup ${alignments} --aligner ${aligner} \\
        --profile ${profile} --min-depth ${params.min_depth} --min-freq ${params.min_freq}
    """

    stub:
    "mkdir pileup"
}

// Reusable by other pipelines (workflows/read-xref): reads (sample, r1, r2), the index and
// the store's members tables (a list, possibly empty).
workflow PILEUP {
    take:
    reads
    index
    store

    main:
    def profiles = QUERY(reads, index).profile
    def refs = REFERENCE(profiles, index, store).ref
    def to_align = refs.join(reads).branch { sample, ref, molecule, r1, r2 ->
        nucleotide: molecule == 'nucleotide'
            return tuple(sample, ref, r1, r2)
        protein: true
            return tuple(sample, ref, r1, r2)
    }
    def alignments = DIAMOND(to_align.protein).alignments
        .map { sample, aligner, r1, r2 -> tuple(sample, aligner, [r1, r2]) }
        .mix(BWAMEM2(to_align.nucleotide).alignments)
    def ref_dirs = refs.map { sample, ref, molecule -> tuple(sample, ref) }
    def pileups = CALL(alignments.join(ref_dirs).join(profiles)).pileup

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
        : INDEX(file(params.targets, checkIfExists: true)).index
    def store = params.store ? channel.fromPath(params.store, checkIfExists: true).collect()
        : channel.value([])
    PILEUP(reads, index, store)
}
