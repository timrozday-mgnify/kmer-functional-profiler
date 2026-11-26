// MGnify proteins -> any HMM collection, by kfpfam's prefilter and HMMER on its candidates,
// by nail's seeds and HMMER (--engine nail), or by HMMER alone (--engine hmmsearch). See README.md. PLAN and FASTA are pfam-reannotate's.

include { PLAN; FASTA; flag } from '../pfam-reannotate/main'
include { NAIL_SEED; NAIL_SEED as NAIL_MEMBERS } from './nail'

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
    if [[ ${'kfpfam' in [params.engine, params.member_filter] || flag(params.profile)} == true ]]; then
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
    publishDir "${params.outdir}/logs", mode: params.publish_mode, pattern: '*.log'

    input:
    tuple val(_arm), val(chunk), path(faa), path(counts), path(seeds)  // seeds: nail's, or none
    path hmm
    path prep
    path index
    path thresholds
    path code

    output:
    tuple val(chunk), path(faa), path("chunk${chunk}.domains.parquet"), emit: domains, optional: true
    path "chunk${chunk}.labels.parquet", emit: labels
    path "chunk${chunk}.stats.tsv", emit: stats
    path "chunk${chunk}.log", emit: log
    tuple val(chunk), path("chunk${chunk}.near.parquet"), emit: near, optional: true

    script:
    """
    set -o pipefail
    ${logged(chunk)}
    ${params.engine == 'kfpfam' ? prefilter(prep, index, faa, chunk, task.cpus) : ''} \\
        ${params.kfpfam_python} annotate.py search --faa ${faa} --pressed models \\
            --thresholds ${thresholds} --threads ${task.cpus} --weak-evalue ${params.weak_evalue} \\
            ${searchable(params.engine, seeds)} ${params.engine in ['hmmsearch', 'release'] ? '' : forward()} \\
            ${flag(params.domains) ? "--domains chunk${chunk}.domains.parquet" : ''} \\
            ${flag(params.members) ? "--near chunk${chunk}.near.parquet --near-slack ${params.member_slack}" : ''} \\
            --labels chunk${chunk}.labels.parquet --stats chunk${chunk}.stats.tsv
    ${searched(chunk, counts)}
    """

    stub:
    "touch chunk${chunk}.domains.parquet chunk${chunk}.labels.parquet chunk${chunk}.log chunk${chunk}.stats.tsv chunk${chunk}.near.parquet"
}

// The filter writes each block's pairs to the pipe as numbers, and HMMER searches them while the
// next block is filtered, as in the raw-reads pipeline's KFPFAM_HMMER. Both use every core.
// PROFILE writes them to a file instead, with each protein's work beside them.
def prefilter(prep, index, faa, chunk, cpus, out = '/dev/stdout', extra = '') {
    def kfpfam = params.kfpfam_modes[params.kfpfam_mode]
    """\
        { time RAYON_NUM_THREADS=${cpus} FILTER_MIN_GAPPED=${kfpfam.min_gapped} \\
        FILTER_GAPPED_FILE=${index}/gapped.txt FILTER_BAND=${params.kfpfam_band} \\
        FILTER_DIAGONALS=${params.kfpfam_diagonals} FILTER_BINARY=1 ${extra} \\
        ${params.kfpfam_filter} join ${prep} ${index}/index.tidx ${kfpfam.bar} ${kfpfam.gate} \\
        ${params.kfpfam_block} ${faa} ${out}; } 2>> chunk${chunk}.log ${out == '/dev/stdout' ? '|' : ''}""".stripIndent()
}

// --engine release: the release's own Pfam table's (protein, family) pairs of the chunk as seeds, so
// HMMER aligns only what that table already has (its domains, with alignments, for lookup-index)
process RELEASE_SEED {
    tag "${chunk}"
    label 'process_small'

    input:
    tuple val(chunk), val(lo), val(hi)
    path thresholds

    output:
    tuple val(chunk), path("chunk${chunk}.seeds.tsv.gz")

    script:
    """
    ${params.python} - <<'EOF'
import duckdb
duckdb.sql(\"\"\"
    COPY (
        SELECT DISTINCT t.name AS model, p.protein_id
        FROM read_parquet('${params.release}/mgy_proteins_pfam.parquet') p
        JOIN read_csv('${thresholds}', delim = '\\t', header = true) t
            ON t.model = format('PF{:05d}', p.pfam_accession)
        WHERE p.protein_id BETWEEN ${lo} AND ${hi - 1} AND p.protein_id % ${params.sample} = 0
    ) TO 'chunk${chunk}.seeds.tsv.gz' (FORMAT csv, DELIMITER '\\t', HEADER, COMPRESSION gzip)
\"\"\")
EOF
    """

    stub:
    "echo | gzip > chunk${chunk}.seeds.tsv.gz"
}

// Where the pairs to search come from, if not stdin
def searchable(engine, seeds) {
    engine in ['nail', 'release'] ?"--seeds ${[seeds].flatten().join(' ')} --block ${params.kfpfam_block}" : engine == 'hmmsearch' ? '--all-pairs' : ''
}

def forward() {
    params.forward_slack == '' ? '' : "--forward-slack ${params.forward_slack}"
}

// The CPU, and the filter's log: its kernel and CPU time, for bench.py. Nextflow's own log has it too.
def logged(chunk) {
    """trap 'cat chunk${chunk}.log >&2' EXIT
    TIMEFORMAT='cpu_s %U %S'
    { grep -m1 'model name' /proc/cpuinfo || uname -m; } > chunk${chunk}.log 2> /dev/null"""
}

// Every protein of the FASTA was read: its counts are the search's
def searched(chunk, counts) {
    """if [[ \$(sed -n 2p chunk${chunk}.stats.tsv | cut -f1,2) != \$(sed -n 2p ${counts}) ]]; then
        echo "chunk${chunk}: proteins or residues searched are not those of ${counts}" >&2
        exit 1
    fi"""
}

// Each protein's cost by hmmsearch and by kfpfam on the same proteins (profile.py)
process PROFILE {
    tag "${chunk}"
    label 'process_search'
    label 'kfpfam'
    publishDir "${params.outdir}/profile", mode: params.publish_mode, pattern: '*.{parquet,json,log}'

    input:
    tuple val(_arm), val(chunk), path(faa), path(_counts)
    path hmm
    path prep
    path index
    path thresholds
    path code

    output:
    path "chunk${chunk}.proteins.parquet", emit: proteins
    path "chunk${chunk}.pairs.parquet", emit: pairs
    path "chunk${chunk}.run.json", emit: run
    path "chunk${chunk}.log", emit: log

    script:
    """
    set -o pipefail
    ${logged(chunk)}
    ${prefilter(prep, index, faa, chunk, task.cpus, "chunk${chunk}.pairs.bin", "FILTER_PER_ORF=chunk${chunk}.orf.tsv")}
    ${params.kfpfam_python} profile.py --faa ${faa} --pressed models --thresholds ${thresholds} \\
        --pairs-bin chunk${chunk}.pairs.bin --per-orf chunk${chunk}.orf.tsv --filter-log chunk${chunk}.log \\
        ${forward()} --threads ${task.cpus} \\
        --proteins chunk${chunk}.proteins.parquet --pairs chunk${chunk}.pairs.parquet --run chunk${chunk}.run.json
    rm chunk${chunk}.pairs.bin chunk${chunk}.orf.tsv
    """

    stub:
    "touch chunk${chunk}.proteins.parquet chunk${chunk}.pairs.parquet chunk${chunk}.run.json chunk${chunk}.log"
}

process PROFILE_SUMMARY {
    label 'process_small'
    publishDir params.outdir, mode: 'copy'

    input:
    path proteins, stageAs: 'proteins/*'
    path pairs, stageAs: 'pairs/*'
    path runs, stageAs: 'runs/*'
    path code

    output:
    path 'profile.*'

    script:
    "${params.python} profile_plot.py summary --proteins proteins --pairs pairs --runs runs --out profile"

    stub:
    "touch profile.bins.tsv"
}

process MEMBER_PAIRS {
    tag "${batch}"
    label 'process_plan'

    input:
    tuple val(batch), val(lo), val(hi), path(near)
    path code

    output:
    path "pairs${batch}", emit: pairs
    path 'shards.txt', emit: shards

    script:
    """
    ${params.python} ${code} --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' pairs --release '${params.release}' \\
        --near ${near} --lo ${lo} --hi ${hi} --sample ${params.sample} --check-sample ${params.check_sample} \\
        ${params.member_filter ? "--all --min-cluster ${params.min_cluster}" : ''} \\
        --width ${params.member_width} --out pairs${batch}
    rm -rf duckdb_tmp
    """

    stub:
    "mkdir pairs${batch}; echo 0 > shards.txt"
}

process MEMBER_FASTA {
    tag "${shard}"
    label 'process_plan'

    input:
    val shard
    path pairs
    path thresholds
    path code

    output:
    tuple val("m${shard}"), path("member${shard}.faa"), path("member${shard}.faa.tsv"), path("member${shard}.pairs")

    script:
    """
    ${params.python} ${code} --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' fasta --release '${params.release}' \\
        --pairs ${pairs} --shard ${shard} --width ${params.member_width} --thresholds ${thresholds} \\
        --block ${params.kfpfam_block} --out member${shard}.faa --pairs-out member${shard}.pairs
    rm -rf duckdb_tmp
    """

    stub:
    "touch member${shard}.faa member${shard}.faa.tsv member${shard}.pairs"
}

process MEMBER_SEARCH {
    tag "${chunk}"
    label 'process_search'
    label 'kfpfam'
    publishDir "${params.outdir}/domains", mode: params.publish_mode, pattern: '*.domains.parquet'
    publishDir "${params.outdir}/labels", mode: params.publish_mode, pattern: '*.labels.parquet'
    publishDir "${params.outdir}/logs", mode: params.publish_mode, pattern: '*.log'

    input:
    tuple val(chunk), path(faa), path(counts), path(pairs), path(seeds)
    path hmm
    path prep
    path index
    path thresholds
    path code

    output:
    tuple val(chunk), path(faa), path("chunk${chunk}.domains.parquet"), emit: domains, optional: true
    path "chunk${chunk}.labels.parquet", emit: labels
    path "chunk${chunk}.stats.tsv", emit: stats
    path "chunk${chunk}.log", emit: log

    script:
    def filtered = params.member_filter as boolean
    """
    set -o pipefail
    ${logged(chunk)}
    ${params.member_filter == 'kfpfam' ? prefilter(prep, index, faa, chunk, task.cpus) : ''} \\
    ${params.kfpfam_python} annotate.py search --faa ${faa} --pressed models \\
        --thresholds ${thresholds} --threads ${task.cpus} --weak-evalue ${params.weak_evalue} \\
        ${filtered ? "--more-pairs ${pairs}" : ''} ${searchable(params.member_filter, seeds)} ${forward()} \\
        ${flag(params.domains) ? "--domains chunk${chunk}.domains.parquet" : ''} \\
        --labels chunk${chunk}.labels.parquet --stats chunk${chunk}.stats.tsv ${filtered ? '' : "< ${pairs}"}
    ${searched(chunk, counts)}
    """

    stub:
    "touch chunk${chunk}.domains.parquet chunk${chunk}.labels.parquet chunk${chunk}.log chunk${chunk}.stats.tsv"
}

process CHECK {
    tag "${chunk}"
    label 'process_search'
    label 'kfpfam'

    input:
    tuple val(name), val(chunk), path(faa), path(domains), val(sample)
    path hmm
    path thresholds
    path code

    output:
    tuple val(name), path("chunk${chunk}.check.tsv")

    script:
    """
    ${params.kfpfam_python} annotate.py check --faa ${faa} --hmm models.h3m \\
        --domains ${domains} --thresholds ${thresholds} --threads ${task.cpus} \\
        --sample ${sample} --check-sample ${params.check_sample} --out chunk${chunk}.check.tsv
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
    if (!(params.engine in ['kfpfam', 'nail', 'hmmsearch', 'release'])) {
        error "--engine ${params.engine}: kfpfam, nail, hmmsearch or release"
    }
    if (params.engine == 'release' && (params.targets != 'all' || flag(params.members))) {
        error "--engine release takes the release table's pairs of every protein: --targets all, no --members"
    }
    if (!(params.kfpfam_mode in params.kfpfam_modes)) {
        error "--kfpfam_mode ${params.kfpfam_mode}: ${params.kfpfam_modes.keySet().join(' or ')}"
    }
    if (flag(params.members) && params.targets != 'reps') {
        error "--members searches representatives first: it needs --targets reps"
    }
    if (!(params.member_filter in ['', 'kfpfam', 'nail'])) {
        error "--member_filter ${params.member_filter}: kfpfam or nail"
    }
    if (params.member_filter && !flag(params.members)) {
        error "--member_filter is for --members true"
    }
    if (params.min_cluster as int > 1 && !params.member_filter) {
        error "--min_cluster leaves the smaller clusters to the prefilter: it needs --member_filter"
    }
    if (!params.hmm) {
        error "--hmm: the collection's HMM file, plain or gzipped"
    }
    if (flag(params.profile) && (params.engine != 'kfpfam' || flag(params.members))) {
        error "--profile times hmmsearch and kfpfam (--kfpfam_mode) on the same proteins: leave --engine kfpfam, without --members"
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
        .ifEmpty { error "nothing to search: no protein is in the sample, or no cluster of --min_cluster" }
    FASTA(ch_chunks)
    if (flag(params.profile)) {
        // a run of its own: both tools on the sample's proteins, timed per protein, and nothing else
        def profiling = [annotate, file("${projectDir}/profile.py", checkIfExists: true)]
        PROFILE(FASTA.out, PREP.out.hmm, PREP.out.prep, PREP.out.index, PREP.out.thresholds, profiling)
        PROFILE_SUMMARY(PROFILE.out.proteins.collect(), PROFILE.out.pairs.collect(), PROFILE.out.run.collect(),
            file("${projectDir}/profile_plot.py", checkIfExists: true))
    } else {
        ch_search = FASTA.out.map { arm, chunk, faa, counts -> [arm, chunk, faa, counts, []] }
        if (params.engine == 'nail') {
            NAIL_SEED(FASTA.out.map { _arm, chunk, faa, _counts -> [chunk, faa] }, PREP.out.hmm)
            ch_search = FASTA.out.map { arm, chunk, faa, counts -> [chunk, arm, faa, counts] }.join(NAIL_SEED.out.seeds)
                .map { chunk, arm, faa, counts, found -> [arm, chunk, faa, counts, found] }
        }
        if (params.engine == 'release') {
            RELEASE_SEED(ch_chunks.map { _arm, chunk, lo, hi, _n -> [chunk, lo, hi] }, PREP.out.thresholds)
            ch_search = FASTA.out.map { arm, chunk, faa, counts -> [chunk, arm, faa, counts] }.join(RELEASE_SEED.out)
                .map { chunk, arm, faa, counts, found -> [arm, chunk, faa, counts, found] }
        }
        SEARCH(ch_search, PREP.out.hmm, PREP.out.prep, PREP.out.index, PREP.out.thresholds, annotate)
        def ch_stats = SEARCH.out.stats
        def ch_check = SEARCH.out.domains.map { chunk, faa, domains -> ['check', chunk, faa, domains, params.sample] }
        if (flag(params.members)) {
            def members = file("${projectDir}/members.py", checkIfExists: true)
            // in the order of the representatives, so that a batch is a range of them, and the batches are every
            // range: the clusters under min_cluster are in no chunk
            ch_near = ch_chunks.map { _arm, chunk, lo, hi, _n -> [chunk, lo, hi] }.join(SEARCH.out.near)
                .toSortedList { a, b -> a[0] <=> b[0] }
                .flatMap { rows ->
                    def parts = rows.groupBy { r -> r[1].intdiv(params.member_batch as long) }.values() as List
                    parts.withIndex().collect { part, i ->
                        [i, i == 0 ? 0L : part[0][1], i + 1 < parts.size() ? parts[i + 1][0][1] : Long.MAX_VALUE,
                         part.collect { r -> r[3] }]
                    }
                }
            MEMBER_PAIRS(ch_near, members)
            ch_shards = MEMBER_PAIRS.out.shards.splitText().map { s -> s.trim() as long }.unique()
            ch_shards.map { s -> "m${s}\t${s * (params.member_width as long)}\t${(s + 1) * (params.member_width as long)}\n" }
                .collectFile(name: 'ranges.members.tsv', seed: "chunk\tlo\thi\n", sort: true, storeDir: params.outdir)
            MEMBER_FASTA(ch_shards, MEMBER_PAIRS.out.pairs.collect(), PREP.out.thresholds, members)
            ch_members = MEMBER_FASTA.out.map { chunk, faa, counts, pairs -> [chunk, faa, counts, pairs, []] }
            if (params.member_filter == 'nail') {
                NAIL_MEMBERS(MEMBER_FASTA.out.map { chunk, faa, _counts, _pairs -> [chunk, faa] }, PREP.out.hmm)
                ch_members = MEMBER_FASTA.out.join(NAIL_MEMBERS.out.seeds)
            }
            MEMBER_SEARCH(ch_members, PREP.out.hmm, PREP.out.prep, PREP.out.index, PREP.out.thresholds, annotate)
            ch_stats = ch_stats.mix(MEMBER_SEARCH.out.stats)
            ch_check = ch_check.mix(MEMBER_SEARCH.out.domains.map { chunk, faa, domains -> ['check.members', chunk, faa, domains, 1] })
        }
        SUMMARY(ch_stats.collect(), tables)
        def ch_checks = channel.empty()
        if (params.check_sample as int > 0) {
            CHECK(ch_check, PREP.out.hmm, PREP.out.thresholds, annotate)
            ch_checks = CHECK.out.groupTuple()
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
}
