// MGnify proteins -> any HMM collection, by kfpfam's prefilter and HMMER on its candidates,
// or by HMMER alone (--engine hmmsearch). See README.md. PLAN and FASTA are pfam-reannotate's.

include { PLAN; FASTA; flag } from '../pfam-reannotate/main'

process PREP {
    label 'process_prep'
    label 'kfpfam'
    publishDir params.outdir, mode: 'copy', pattern: '{thresholds.tsv,hmm.sha256}'

    input:
    path hmm
    path table, stageAs: 'table/*'
    path code

    output:
    path 'models.{hmm,h3m,h3i,h3f,h3p}', emit: hmm
    path 'prep', emit: prep
    path 'index', emit: index  // the filter's tiled index and gapped bars; empty for hmmsearch
    path 'thresholds.tsv', emit: thresholds
    path 'hmm.sha256', emit: sha256

    script:
    def kfpfam = params.kfpfam_modes[params.kfpfam_mode]
    def source = table ? "--table ${table}" : "--cutoff ${params.cutoff}"
    """
    if [[ ${hmm} == *.gz ]]; then gunzip -c ${hmm} > tmp.hmm; else cp -L ${hmm} tmp.hmm; fi
    mv tmp.hmm models.hmm
    ${params.kfpfam_python} ${params.kfpfam_prep_py} models.hmm - prep
    ${params.kfpfam_python} tables.py thresholds --hmm models.hmm --models prep/models.tsv ${source}
    ${params.kfpfam_python} annotate.py press --hmm models.hmm --models prep/models.tsv --out models
    mkdir index
    if [[ ${params.engine} == kfpfam ]]; then
        export RAYON_NUM_THREADS=${task.cpus}
        ${params.kfpfam_python} seeds.py bars --models prep/models.tsv --thresholds thresholds.tsv \\
            --ga ga.tsv --gapped index/gapped.txt --slack ${kfpfam.slack}
        ${params.kfpfam_filter} seeds prep ga.tsv seeds.f32 ${params.kfpfam_seed}
        ${params.kfpfam_python} seeds.py choose --table seeds.f32 --shape ${params.kfpfam_seed} \\
            --theta ${kfpfam.theta} --out T.txt
        FILTER_T_FILE=T.txt ${params.kfpfam_filter} build prep ${params.kfpfam_seed} 0 index.idx
        ${params.kfpfam_filter} tile prep index.idx index/index.tidx
        rm index.idx seeds.f32
    fi
    """

    stub:
    "mkdir prep index; touch models.hmm models.h3m models.h3i models.h3f models.h3p thresholds.tsv hmm.sha256"
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
    def kfpfam = params.kfpfam_modes[params.kfpfam_mode]
    def prefilter = """\
        RAYON_NUM_THREADS=${task.cpus} FILTER_MIN_GAPPED=${kfpfam.min_gapped} \\
        FILTER_GAPPED_FILE=${index}/gapped.txt FILTER_BAND=${params.kfpfam_band} \\
        FILTER_DIAGONALS=${params.kfpfam_diagonals} FILTER_BINARY=1 \\
        ${params.kfpfam_filter} join ${prep} ${index}/index.tidx ${kfpfam.bar} ${kfpfam.gate} \\
        ${params.kfpfam_block} ${faa} /dev/stdout |""".stripIndent()
    """
    set -o pipefail
    ${params.engine == 'kfpfam' ? prefilter : ''} \\
        ${params.kfpfam_python} annotate.py search --faa ${faa} --pressed models \\
            --thresholds ${thresholds} --threads ${task.cpus} --weak-evalue ${params.weak_evalue} \\
            ${params.engine == 'kfpfam' ? '' : '--all-pairs'} \\
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

process COMPARE {
    tag "${chunk}"
    label 'process_small'

    input:
    tuple val(chunk), val(lo), val(hi), path(domains)
    path code

    output:
    path "chunk${chunk}.compare.tsv"

    script:
    """
    ${params.python} ${code} --domains ${domains} --table '${params.release}/mgy_proteins_pfam.parquet' \\
        --lo ${lo} --hi ${hi} --sample ${params.sample} --threads ${task.cpus} --out chunk${chunk}.compare.tsv
    """

    stub:
    "touch chunk${chunk}.compare.tsv"
}

process SUMMARY {
    label 'process_small'
    publishDir params.outdir, mode: 'copy'

    input:
    path stats, stageAs: 'chunks/*'
    path code

    output:
    path 'summary.tsv'

    script:
    "${params.python} tables.py summary chunks"

    stub:
    "touch summary.tsv"
}

process CHECK_SUMMARY {
    tag "${name}"
    label 'process_small'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(name), path(checks, stageAs: 'chunks/*')
    path code

    output:
    path "${name}.tsv"

    // by directory: thousands of names would not fit on a command line
    script:
    "${params.python} tables.py check-summary chunks --out ${name}.tsv"

    stub:
    "touch ${name}.tsv"
}

workflow {
    if (!(params.targets in ['reps', 'all'])) {
        error "--targets ${params.targets}: reps or all"
    }
    if (!(params.engine in ['kfpfam', 'hmmsearch'])) {
        error "--engine ${params.engine}: kfpfam or hmmsearch"
    }
    if (!(params.kfpfam_mode in params.kfpfam_modes)) {
        error "--kfpfam_mode ${params.kfpfam_mode}: ${params.kfpfam_modes.keySet().join(' or ')}"
    }
    if (!params.hmm) {
        error "--hmm: the collection's HMM file, plain or gzipped"
    }
    if (params.check_sample as int > 0 && !flag(params.domains)) {
        error "--check_sample compares domains: it needs --domains true"
    }
    if (flag(params.compare_release) && !(flag(params.domains) && params.targets == 'all')) {
        error "--compare_release compares domains of every protein: it needs --domains true --targets all"
    }
    if (workflow.containerEngine && !params.kfpfam_sif) {
        error "--kfpfam_sif: kfpfam's image (kfpfam's docker/filter.Dockerfile)"
    }
    def table = params.thresholds ? file(params.thresholds, checkIfExists: true) : []
    // inputs, so that containers see them and -resume reruns what an edit changes
    def annotate = file("${projectDir}/annotate.py", checkIfExists: true)
    def tables = file("${projectDir}/tables.py", checkIfExists: true)
    def seeds = file("${projectDir}/seeds.py", checkIfExists: true)
    PREP(file(params.hmm, checkIfExists: true), table, [annotate, tables, seeds])
    ch_chunks = PLAN().splitCsv(sep: '\t', header: true)
        .map { r -> ['real', r.chunk as int, r.lo as long, r.hi as long, r.proteins as long] }
    FASTA(ch_chunks)
    SEARCH(FASTA.out, PREP.out.hmm, PREP.out.prep, PREP.out.index, PREP.out.thresholds, annotate)
    SUMMARY(SEARCH.out.stats.collect(), tables)
    def ch_checks = channel.empty()
    if (params.check_sample as int > 0) {
        CHECK(SEARCH.out.domains, PREP.out.hmm, PREP.out.thresholds, annotate)
        ch_checks = CHECK.out.collect().map { c -> ['check', c] }
    }
    if (flag(params.compare_release)) {
        def compare = file("${projectDir}/compare.py", checkIfExists: true)
        COMPARE(
            ch_chunks.map { _arm, chunk, lo, hi, _n -> [chunk, lo, hi] }
                .join(SEARCH.out.domains.map { chunk, _faa, domains -> [chunk, domains] }),
            compare,
        )
        ch_checks = ch_checks.mix(COMPARE.out.collect().map { c -> ['compare', c] })
    }
    CHECK_SUMMARY(ch_checks, tables)
}
