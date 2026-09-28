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
    publishDir params.outdir, mode: 'copy', pattern: 'index_*/meta.json'

    input:
    tuple val(name), val(args)
    path members

    output:
    tuple val(name), path("index_${name}"), emit: index
    path "index_${name}/meta.json"

    script:
    "${params.kfp} index ${members} index_${name} ${args}"

    stub:
    "mkdir index_${name} && touch index_${name}/meta.json"
}

process IMPORT_SKETCHES {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', pattern: 'index_*/meta.json'

    input:
    path sketches

    output:
    tuple val('fmh_compat'), path('index_fmh_compat'), emit: index
    path 'index_fmh_compat/meta.json'

    script:
    "${params.kfp} import-sourmash ${sketches} index_fmh_compat --ksize ${params.ksize}"

    stub:
    "mkdir index_fmh_compat && touch index_fmh_compat/meta.json"
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
    """
    ${params.venv}/bin/iss generate --genomes ${fna} --model ${params.iss_model} \\
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
    publishDir "${params.outdir}/profiles", mode: 'copy', saveAs: { "seed${seed}_${name}.tsv" }

    input:
    tuple val(seed), path(r1), path(r2), val(name), path(index)

    output:
    tuple val(seed), val(name), path('profile.tsv'), emit: profile

    script:
    "${params.kfp} query ${index} ${r1} ${r2} --out profile.tsv"

    stub:
    "touch profile.tsv"
}

process SCORE {
    tag "seed ${seed} ${name}"
    label 'process_single'

    input:
    tuple val(seed), val(name), path(profile), path(truth)

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

workflow {
    FETCH()
    MEMBERS(FETCH.out.faa, FETCH.out.kos)
    ch_configs = channel.fromList(params.indexes).map { cfg -> [cfg.name, cfg.args] }
    INDEX(ch_configs, MEMBERS.out.members)
    IMPORT_SKETCHES(FETCH.out.sketches)
    ch_indexes = INDEX.out.index.mix(IMPORT_SKETCHES.out.index)

    SAMPLE(channel.of(1..params.replicates), FETCH.out.genomes)
    SIMULATE(SAMPLE.out.sample.map { seed, fna, _genes -> [seed, fna] })
    TRUTH(SAMPLE.out.sample.join(SIMULATE.out.reads), FETCH.out.kos)

    PROFILE(SIMULATE.out.reads.combine(ch_indexes))
    SCORE(PROFILE.out.profile.combine(TRUTH.out.truth, by: 0))
    SUMMARY(SCORE.out.score.collect())
}
