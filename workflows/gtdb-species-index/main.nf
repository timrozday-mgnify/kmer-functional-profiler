// GTDB species index (plan: phase 11, step 11; README.md): a GTDB release's genomes, at most
// --max_per_species per species (representative first), their DNA from NCBI, genes called by
// pyrodigal, annotated with the index (annotate-genomes) and aggregated by GTDB species
// (species-index --genomes).

process METADATA {
    label 'process_single'
    storeDir params.data_dir

    output:
    path 'metadata/*', emit: metadata

    script:
    """
    mkdir metadata
    for f in ${params.metadata.join(' ')}; do
        curl -fsSL --retry 3 -o metadata/\$f ${params.gtdb_url}/\$f
    done
    """

    stub:
    "mkdir metadata && touch metadata/${params.metadata[0]}"
}

process PICK {
    label 'process_single'
    publishDir params.outdir, mode: 'copy'

    input:
    path metadata
    path code, stageAs: 'code/*'  // gtdb.py: only here so -resume reruns on changes

    output:
    path 'genomes.tsv', emit: genomes
    path 'annotate.tsv', emit: annotate

    script:
    """
    ${params.python} ${code} ${metadata} --max-per-species ${params.max_per_species} \\
        --min-completeness ${params.min_completeness} --max-contamination ${params.max_contamination} \\
        ${params.species ? "--species ${params.species}" : ''} \\
        ${params.max_species ? "--max-species ${params.max_species}" : ''} \\
        --seed ${params.seed} --ncbi ${params.ncbi_url}
    """

    stub:
    "printf 'genome\\turl\\n' > genomes.tsv && printf 'genome\\tpath\\n' > annotate.tsv"
}

process FETCH {
    tag "${shard}"
    label 'process_single'
    maxForks params.max_downloads

    input:
    tuple val(shard), val(rows)

    output:
    tuple val(shard), path('dna/*'), emit: dna

    script:
    // one URL per genome; a genome NCBI no longer serves fails the shard (rerun with it
    // dropped from the metadata)
    def list = rows.collect { r -> "${r.genome}\t${r.url}" }.join('\n')
    """
    mkdir dna
    printf '%s\\n' '${list}' | while IFS=\$'\\t' read -r genome url; do
        curl -fsSL --retry 3 -o dna/\${genome}.fna.gz "\${url}"
    done
    """

    stub:
    "mkdir dna && touch dna/${rows[0].genome}.fna.gz"
}

process GENES {
    tag "${shard}"
    label 'process_medium'
    container 'quay.io/biocontainers/pyrodigal:3.7.1--py312h247cb63_1'
    publishDir "${params.outdir}", mode: 'copy', enabled: params.publish_proteins

    input:
    tuple val(shard), path(dna, stageAs: 'dna/*')

    output:
    path 'faa/*', emit: faa

    script:
    // single mode trains on each genome; genomes too short to train on (< 20 kb) fall back to
    // meta. Table 11 throughout; Prodigal's stop '*' dropped.
    """
    mkdir faa
    for f in dna/*.fna.gz; do
        g=\$(basename "\$f" .fna.gz)
        gzip -dc "\$f" > genome.fna
        pyrodigal -p single -i genome.fna -a genes.faa -o /dev/null 2>/dev/null \\
            || pyrodigal -p meta -i genome.fna -a genes.faa -o /dev/null
        sed 's/\\*\$//' genes.faa > faa/\${g}.faa
    done
    rm -f genome.fna genes.faa
    """

    stub:
    "mkdir faa && touch faa/x.faa"
}

process ANNOTATE {
    label 'process_high_memory'

    input:
    path index
    path annotate
    path faa, stageAs: 'faa/*'
    path code, stageAs: 'code*/*'

    output:
    path 'genome_index', emit: index

    script:
    "${params.kfp} annotate-genomes ${index} ${annotate} genome_index"

    stub:
    "mkdir genome_index"
}

process SPECIES_INDEX {
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy'

    input:
    path index
    path genome_index
    path code, stageAs: 'code*/*'

    output:
    path 'species_index'
    path 'species_index.json'

    script:
    "${params.kfp} species-index ${index} species_index --genomes ${genome_index} ${params.species_index_args} > species_index.json"

    stub:
    "mkdir species_index && touch species_index.json"
}

process PANEL {
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy'

    input:
    path index
    path genome_index
    path code, stageAs: 'code*/*'

    output:
    path 'panel'
    path 'panel.json'

    script:
    "${params.kfp_genomes} panel ${index} panel --genomes ${genome_index} ${params.panel_args} > panel.json"

    stub:
    "mkdir panel && touch panel.json"
}

workflow {
    if (!params.index) {
        error "--index: the index to annotate with (a profile's index)"
    }
    def ch_index = channel.value(file(params.index, checkIfExists: true))
    // ponytail: tracks the Python package only; Rust kernel changes still need a fresh run
    def ch_code = channel.fromPath("${projectDir}/../../python/{kmer_functional_profiler,kfp_genomes}/*.py").collect()
    PICK(METADATA().metadata, file("${projectDir}/gtdb.py"))
    def ch_shards = PICK.out.genomes
        .splitCsv(header: true, sep: '\t')
        .buffer(size: params.shard_size, remainder: true)
        .map { rows -> [rows[0].genome, rows] }  // a shard is named by its first genome
    GENES(FETCH(ch_shards).dna)
    ANNOTATE(ch_index, PICK.out.annotate, GENES.out.faa.collect(), ch_code)
    SPECIES_INDEX(ch_index, ANNOTATE.out.index, ch_code)
    if (params.panel) {
        PANEL(ch_index, ANNOTATE.out.index, ch_code)  // the genome panel (phase 12, step 2)
    }
}
