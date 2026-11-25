// MGnify proteins -> any HMM collection, by kfpfam's prefilter and HMMER on its candidates.
// See README.md. PLAN and FASTA are pfam-reannotate's.

include { PLAN; FASTA; flag } from '../pfam-reannotate/main'

process PREP {
    label 'process_prep'
    label 'kfpfam'
    publishDir params.outdir, mode: 'copy', pattern: 'thresholds.tsv'

    input:
    path hmm
    path table, stageAs: 'table/*'
    path code

    output:
    path 'models.{hmm,h3m,h3i,h3f,h3p}', emit: hmm
    path 'prep', emit: prep
    path 'index.tidx', emit: index
    path 'thresholds.tsv', emit: thresholds

    script:
    def source = table ? "--table ${table}" : "--cutoff ${params.cutoff}"
    """
    if [[ ${hmm} == *.gz ]]; then gunzip -c ${hmm} > tmp.hmm; else cp -L ${hmm} tmp.hmm; fi
    mv tmp.hmm models.hmm
    ${params.kfpfam_python} ${params.kfpfam_prep_py} models.hmm - prep
    ${params.kfpfam_python} tables.py thresholds --hmm models.hmm --models prep/models.tsv ${source}
    ${params.kfpfam_python} annotate.py press --hmm models.hmm --models prep/models.tsv --out models
    ${params.kfpfam_filter} build prep 4 ${params.kfpfam_bits} index.idx
    ${params.kfpfam_filter} tile prep index.idx index.tidx
    rm index.idx
    """

    stub:
    "mkdir prep; touch models.hmm models.h3m models.h3i models.h3f models.h3p index.tidx thresholds.tsv"
}

process SEARCH {
    tag "${chunk}"
    label 'process_search'
    label 'kfpfam'
    publishDir "${params.outdir}/domains", mode: params.publish_mode, pattern: '*.domains.parquet'
    publishDir "${params.outdir}/labels", mode: params.publish_mode, pattern: '*.labels.parquet'

    input:
    tuple val(arm), val(chunk), path(faa), path(counts)
    path hmm
    path prep
    path index
    path thresholds
    path code

    output:
    tuple val(chunk), path(faa), path("chunk${chunk}.domains.parquet"), emit: domains, optional: true
    path "chunk${chunk}.labels.parquet", emit: labels
    path "chunk${chunk}.stats.tsv", emit: stats

    // The filter writes each block's pairs to the pipe as numbers, and HMMER searches them while the
    // next block is filtered, as in the raw-reads pipeline's KFPFAM_HMMER. Both use every core.
    script:
    """
    set -o pipefail
    RAYON_NUM_THREADS=${task.cpus} FILTER_MIN_GAPPED=${params.kfpfam_min_gapped} FILTER_BINARY=1 \\
        ${params.kfpfam_filter} join ${prep} ${index} ${params.kfpfam_bar} ${params.kfpfam_gate} \\
        ${params.kfpfam_block} ${faa} /dev/stdout \\
        | ${params.kfpfam_python} annotate.py search --faa ${faa} --pressed models \\
            --thresholds ${thresholds} --threads ${task.cpus} --weak-evalue ${params.weak_evalue} \\
            ${flag(params.domains) ? "--domains chunk${chunk}.domains.parquet" : ''} \\
            --labels chunk${chunk}.labels.parquet --stats chunk${chunk}.stats.tsv
    if [[ \$(sed -n 2p chunk${chunk}.stats.tsv | cut -f1,2) != \$(sed -n 2p ${counts}) ]]; then
        echo "chunk${chunk}: proteins or residues searched are not those of ${counts}" >&2
        exit 1
    fi
    """

    stub:
    "touch chunk${chunk}.domains.parquet chunk${chunk}.labels.parquet chunk${chunk}.stats.tsv"
}

process CHECK {
    tag "${chunk}"
    label 'process_search'
    label 'kfpfam'

    input:
    tuple val(chunk), path(faa), path(domains)
    path hmm
    path thresholds
    path code

    output:
    path "chunk${chunk}.check.tsv"

    script:
    """
    ${params.kfpfam_python} annotate.py check --faa ${faa} --hmm models.h3m \\
        --domains ${domains} --thresholds ${thresholds} --threads ${task.cpus} \\
        --sample ${params.sample} --check-sample ${params.check_sample} --out chunk${chunk}.check.tsv
    """

    stub:
    "touch chunk${chunk}.check.tsv"
}

process SUMMARY {
    label 'process_small'
    publishDir params.outdir, mode: 'copy'

    input:
    path stats
    path code

    output:
    path 'summary.tsv'

    script:
    "${params.python} tables.py summary ${stats}"

    stub:
    "touch summary.tsv"
}

process CHECK_SUMMARY {
    label 'process_small'
    publishDir params.outdir, mode: 'copy'

    input:
    path checks
    path code

    output:
    path 'check.tsv'

    script:
    "${params.python} tables.py check-summary ${checks}"

    stub:
    "touch check.tsv"
}

workflow {
    if (!(params.targets in ['reps', 'all'])) {
        error "--targets ${params.targets}: reps or all"
    }
    if (!params.hmm) {
        error "--hmm: the collection's HMM file, plain or gzipped"
    }
    if (params.check_sample as int > 0 && !flag(params.domains)) {
        error "--check_sample compares domains: it needs --domains true"
    }
    if (workflow.containerEngine && !params.kfpfam_sif) {
        error "--kfpfam_sif: kfpfam's image (kfpfam's docker/filter.Dockerfile)"
    }
    def table = params.thresholds ? file(params.thresholds, checkIfExists: true) : []
    // inputs, so that containers see them and -resume reruns what an edit changes
    def annotate = file("${projectDir}/annotate.py", checkIfExists: true)
    def tables = file("${projectDir}/tables.py", checkIfExists: true)
    PREP(file(params.hmm, checkIfExists: true), table, [annotate, tables])
    ch_chunks = PLAN().splitCsv(sep: '\t', header: true)
        .map { r -> ['real', r.chunk as int, r.lo as long, r.hi as long, r.proteins as long] }
    FASTA(ch_chunks)
    SEARCH(FASTA.out, PREP.out.hmm, PREP.out.prep, PREP.out.index, PREP.out.thresholds, annotate)
    SUMMARY(SEARCH.out.stats.collect(), tables)
    if (params.check_sample as int > 0) {
        CHECK(SEARCH.out.domains, PREP.out.hmm, PREP.out.thresholds, annotate)
        CHECK_SUMMARY(CHECK.out.collect(), tables)
    }
}
