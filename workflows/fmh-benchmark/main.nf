// fmh-funprofiler benchmark: KO detection by kmer-functional-profiler indexes and by the
// fmh-funprofiler KO sketches (compat mode) on InSilicoSeq metagenomes. See README.md.


process FETCH {
    label 'process_single'
    storeDir params.data_dir

    output:
    path 'genomes_extracted_from_kegg', emit: genomes
    path 'protein_ref_db_giant.faa',    emit: faa
    path 'present_genes_and_koids.csv', emit: kos
    path "${params.sketches}",          emit: sketches

    script:
    """
    Z=https://zenodo.org/records/10055954/files
    curl -fsSL -o genomes.zip "\$Z/genomes_extracted_from_kegg.zip?download=1"
    unzip -q genomes.zip -x '*.DS_Store' '__MACOSX/*' && rm genomes.zip
    for f in protein_ref_db_giant.faa present_genes_and_koids.csv ${params.sketches}; do
        curl -fsSL -o "\$f" "\$Z/\$f?download=1"
    done
    """

    stub:
    "mkdir genomes_extracted_from_kegg && touch protein_ref_db_giant.faa present_genes_and_koids.csv ${params.sketches}"
}

process MEMBERS {
    label 'process_medium'

    input:
    path faa
    path kos

    output:
    path 'members.parquet', emit: members

    script:
    "${params.bench} members --faa ${faa} --kos ${kos}"

    stub:
    "touch members.parquet"
}

process INDEX {
    tag "${name}"
    label 'process_high_memory'
    // meta.json is copied out of the index: a file inside a declared directory output is
    // folded into that directory and never matches the publish pattern
    publishDir params.outdir, mode: 'copy', pattern: 'meta.json', saveAs: { "index_${name}/meta.json" }

    input:
    tuple val(name), val(args)
    path members
    path code, stageAs: 'code/*'  // package sources: only here so -resume rebuilds on changes

    output:
    tuple val(name), path("index_${name}"), emit: index
    path 'meta.json'

    script:
    "${params.kfp} index ${members} index_${name} ${args} && cp index_${name}/meta.json meta.json"

    stub:
    "mkdir index_${name} && touch index_${name}/meta.json meta.json"
}

process IMPORT_SKETCHES {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', pattern: 'meta.json', saveAs: { 'index_fmh_compat/meta.json' }

    input:
    path sketches
    path code, stageAs: 'code/*'  // package sources: only here so -resume rebuilds on changes

    output:
    tuple val('fmh_compat'), path('index_fmh_compat'), emit: index
    path 'meta.json'

    script:
    """
    ${params.kfp} import-sourmash ${sketches} index_fmh_compat --ksize ${params.ksize}
    cp index_fmh_compat/meta.json meta.json
    """

    stub:
    "mkdir index_fmh_compat && touch index_fmh_compat/meta.json meta.json"
}

process SAMPLE {
    tag "seed ${seed}"
    label 'process_single'

    input:
    val seed
    path genomes

    output:
    tuple val(seed), path('sample.fna'), path('genes.parquet'), emit: sample
    path 'genomes.txt'

    script:
    "${params.bench} sample --genomes-dir ${genomes} --n ${params.n_genomes} --seed ${seed}"

    stub:
    "touch sample.fna genes.parquet genomes.txt"
}

process SIMULATE {
    tag "seed ${seed}"
    label 'process_medium'

    input:
    tuple val(seed), path(fna)

    output:
    tuple val(seed), path('reads_R1.fastq.gz'), path('reads_R2.fastq.gz'), emit: reads

    script:
    // bench.py iss patches the perfect model's bugs in iss 2.0.1
    """
    ${params.iss_mode == 'perfect' ? "${params.bench} iss generate --mode perfect" : "${params.venv}/bin/iss generate --model ${params.iss_model}"} --genomes ${fna} \\
        --n_reads ${params.n_reads} --abundance lognormal --seed ${seed} \\
        --cpus ${task.cpus} --compress --output reads
    """

    stub:
    "touch reads_R1.fastq.gz reads_R2.fastq.gz"
}

process TRUTH {
    tag "seed ${seed}"
    label 'process_medium'
    publishDir "${params.outdir}/truth", mode: 'copy', saveAs: { "seed${seed}.csv" }

    input:
    tuple val(seed), path(fna), path(genes), path(r1), path(r2)
    path kos
    path code, stageAs: 'code/*'  // bench.py: only here so -resume reruns on changes

    output:
    tuple val(seed), path('truth.csv'), emit: truth

    script:
    """
    ${params.bench} truth --fna ${fna} --genes ${genes} --kos ${kos} --r1 ${r1} --r2 ${r2} \\
        --threads ${task.cpus}
    """

    stub:
    "touch truth.csv"
}

process PROFILE {
    tag "seed ${seed} ${name}"
    label 'process_single'
    // one saveAs for both outputs: without the per-file branch they overwrite each other
    publishDir params.outdir, mode: 'copy', saveAs: { f ->
        f.endsWith('.parquet') ? "kmers/seed${seed}_${name}.parquet" : "profiles/seed${seed}_${name}.tsv"
    }

    input:
    tuple val(seed), path(r1), path(r2), val(name), path(index)
    path code, stageAs: 'code/*'  // query sources: only here so -resume reruns on changes

    output:
    tuple val(seed), val(name), path('profile.tsv'), path('kmers.parquet'), emit: profile

    script:
    "${params.kfp} query ${index} ${r1} ${r2} --out profile.tsv --draws ${params.draws} --kmers kmers.parquet"

    stub:
    "touch profile.tsv kmers.parquet"
}

process SCORE {
    tag "seed ${seed} ${name}"
    label 'process_single'

    input:
    tuple val(seed), val(name), path(profile), path(truth)
    path code, stageAs: 'code/*'  // bench.py: only here so -resume reruns on changes

    output:
    path 'score.tsv', emit: score

    script:
    """
    ${params.bench} score --truth ${truth} --profile ${profile} --sample seed${seed} --index ${name} \\
        --min-hits ${params.min_hits.toString().tokenize(',').join(' ')}
    """

    stub:
    "touch score.tsv"
}

process DETECTED {
    tag "seed ${seed} ${name}"
    label 'process_single'
    publishDir "${params.outdir}/detected", mode: 'copy', saveAs: { "seed${seed}_${name}.tsv" }

    input:
    tuple val(name), val(seed), path(profile), path(kmers), path(truth), path(fna), path(index)
    path code, stageAs: 'code/*'  // bench.py: only here so -resume reruns on changes

    output:
    path 'detected.tsv', emit: detected

    script:
    """
    ${params.bench} detected --truth ${truth} --profile ${profile} --kmers ${kmers} \\
        --index-dir ${index} --fna ${fna} --sample seed${seed} --index ${name}
    """

    stub:
    "touch detected.tsv"
}

process SUMMARY {
    label 'process_single'
    publishDir params.outdir, mode: 'copy'

    input:
    path scores, stageAs: 'score*.tsv'

    output:
    path 'summary.tsv'
    path 'scores.tsv'

    script:
    "${params.bench} summary ${scores}"

    stub:
    "touch summary.tsv scores.tsv"
}

// Output of a git command in the pipeline's checkout (for run.json).
def git(List cmd) {
    return (['git', '-C', projectDir.toString()] + cmd).execute().text.trim()
}

workflow {
    main:
    FETCH()
    MEMBERS(FETCH.out.faa, FETCH.out.kos)
    // ponytail: tracks the Python package only; Rust kernel changes still need a fresh run
    ch_code = channel.fromPath("${projectDir}/../../python/kmer_functional_profiler/*.py").collect()
    ch_configs = channel.fromList(params.indexes).map { cfg -> [cfg.name, cfg.args] }
    INDEX(ch_configs, MEMBERS.out.members, ch_code)
    IMPORT_SKETCHES(FETCH.out.sketches, ch_code)
    ch_indexes = INDEX.out.index.mix(IMPORT_SKETCHES.out.index)

    SAMPLE(channel.of(1..params.replicates), FETCH.out.genomes)
    SIMULATE(SAMPLE.out.sample.map { seed, fna, _genes -> [seed, fna] })
    TRUTH(SAMPLE.out.sample.join(SIMULATE.out.reads), FETCH.out.kos, file("${projectDir}/bench.py"))

    PROFILE(SIMULATE.out.reads.combine(ch_indexes), ch_code)
    ch_scored = PROFILE.out.profile.combine(TRUTH.out.truth, by: 0)  // seed, name, profile, kmers, truth
    SCORE(ch_scored.map { seed, name, profile, _kmers, truth -> [seed, name, profile, truth] }, file("${projectDir}/bench.py"))
    DETECTED(
        ch_scored
            .combine(SAMPLE.out.sample.map { seed, fna, _genes -> [seed, fna] }, by: 0)
            .map { seed, name, profile, kmers, truth, fna -> [name, seed, profile, kmers, truth, fna] }
            .combine(ch_indexes, by: 0),
        file("${projectDir}/bench.py"),
    )
    DETECTED.out.detected.collectFile(name: 'detected.tsv', keepHeader: true, storeDir: params.outdir)
    SUMMARY(SCORE.out.score.collect())

    onComplete:
    // run.json: what produced the results in outdir (reads, indexes, code, status)
    def reads = params.iss_mode == 'perfect' ? 'error-free (iss perfect)' : "iss ${params.iss_model} (${params.iss_mode})"
    def info = [
        description: "${params.replicates} metagenomes x ${params.n_genomes} KEGG genomes, " +
            "${params.n_reads} reads, ${reads}; ${params.indexes.size()} indexes + fmh_compat; draws ${params.draws}",
        success: workflow.success,
        exit_status: workflow.exitStatus,
        start: workflow.start.toString(),
        complete: workflow.complete.toString(),
        duration: workflow.duration.toString(),
        command_line: workflow.commandLine,
        profile: workflow.profile,
        resume: workflow.resume,
        session_id: workflow.sessionId.toString(),
        run_name: workflow.runName,
        nextflow: workflow.nextflow.version.toString(),
        code: [commit: git(['rev-parse', 'HEAD']), branch: git(['rev-parse', '--abbrev-ref', 'HEAD']),
               uncommitted_changes: !git(['status', '--porcelain']).isEmpty()],
        params: params,
    ]
    def out = file(params.outdir)
    out.mkdirs()
    out.resolve('run.json').text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info)) + '\n'
}
