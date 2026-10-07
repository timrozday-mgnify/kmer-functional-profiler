// Species model strain hold-out benchmark (plan: phase 11, step 10; README.md): samples of
// non-representative MGnify catalogue genomes, simulated with InSilicoSeq and profiled; the
// species index is built with those genomes left out, so every sample strain is new to it.
// Arms: the species model and its ablations, genome mode (G1) on the representatives with
// kfp-prior (B1) on top, and sylph on the representatives' DNA.

// --name FILE for an input given and not empty (an arm without that output passes [])
def opt(name, files) {
    def given = (files instanceof List ? files : [files]).findAll { f -> f.size() > 0 }
    return given ? "--${name} ${given[0]}" : ''
}

process METADATA {
    label 'process_single'
    storeDir params.data_dir

    output:
    path 'genomes-all_metadata.tsv', emit: metadata

    script:
    "curl -fsSL --retry 3 -o genomes-all_metadata.tsv ${params.catalogue_url}/genomes-all_metadata.tsv"

    stub:
    "touch genomes-all_metadata.tsv"
}

process PICK {
    label 'process_single'
    publishDir params.outdir, mode: 'copy'

    input:
    path metadata
    path code, stageAs: 'code*/*'  // bench.py: only here so -resume reruns on changes

    output:
    path 'samples.tsv', emit: samples
    path 'exclude.txt', emit: exclude
    path 'species.txt', emit: species

    script:
    """
    ${params.bench} uhgg-pick --metadata ${metadata} --replicates ${params.replicates} \\
        --per-sample ${params.per_sample} --min-genomes ${params.min_genomes} \\
        --min-completeness ${params.min_completeness} --max-contamination ${params.max_contamination} \\
        --distractors ${params.distractors} --depth-min ${params.depth_min} --depth-max ${params.depth_max} \\
        --seed ${params.seed}
    """

    stub:
    "printf 'sample\\tgenome\\tspecies\\tdepth\\n' > samples.tsv && touch exclude.txt species.txt"
}

process FETCH_CATALOGUE {
    label 'process_medium'

    input:
    path metadata
    path species

    output:
    path 'catalogue', emit: catalogue

    script:
    """
    ${projectDir}/fetch.sh species ${params.catalogue_url} ${species} catalogue ${task.cpus}
    cp ${metadata} catalogue/genomes-all_metadata.tsv
    """

    stub:
    "mkdir catalogue"
}

process FETCH_GENOMES {
    label 'process_medium'

    input:
    path samples

    output:
    path 'genomes', emit: genomes

    script:
    "${projectDir}/fetch.sh genomes ${params.catalogue_url} ${samples} genomes ${task.cpus}"

    stub:
    "mkdir genomes"
}

process SIMULATE {
    tag "sample ${sample}"
    label 'process_medium'

    input:
    val sample
    path samples
    path genomes
    path code, stageAs: 'code*/*'

    output:
    tuple val(sample), path('reads_R1.fastq.gz'), path('reads_R2.fastq.gz'), emit: reads

    script:
    // bench.py iss patches the perfect model's bugs in iss 2.0.1
    def iss = params.iss_mode == 'perfect' ? "${params.bench} iss generate --mode perfect" : "${params.venv}/bin/iss generate --model ${params.iss_model}"
    """
    ${params.bench} uhgg-sample --samples ${samples} --sample ${sample} --genomes-dir ${genomes}
    # iss 2.0.1 splits --coverage_file work by its default --n_reads, so cpus whose chunk
    # gets no reads write no temporary file and the final concatenation fails: create them
    for i in \$(seq 0 ${task.cpus - 1}); do touch reads.iss.tmp.\${i}_R1.fastq reads.iss.tmp.\${i}_R2.fastq; done
    ${iss} --genomes sample.fna --coverage_file coverage.txt --seed ${sample} \\
        --cpus ${task.cpus} --compress --output reads
    """

    stub:
    "touch reads_R1.fastq.gz reads_R2.fastq.gz"
}

process INDEX {
    label 'process_high_memory'

    input:
    path members, stageAs: 'members/*'
    path pfam, stageAs: 'pfam/*'
    path code, stageAs: 'code*/*'

    output:
    path 'index', emit: index

    script:
    "${params.kfp} index ${members} index ${pfam ? "--pfam ${pfam}" : ''} ${params.index_args}"

    stub:
    "mkdir index"
}

process PROFILE {
    tag "sample ${sample}"
    label 'process_medium'
    publishDir "${params.outdir}/profiles", mode: 'copy', saveAs: { "sample${sample}.tsv" }

    input:
    tuple val(sample), path(r1), path(r2)
    path index
    path code, stageAs: 'code*/*'

    output:
    tuple val(sample), path('profile.tsv'), emit: profile

    script:
    "${params.kfp} query ${index} ${r1} ${r2} --out profile.tsv ${params.query_args}"

    stub:
    "touch profile.tsv"
}

process SPECIES_INDEX {
    tag "${name}"
    label 'process_high_memory'

    input:
    tuple val(name), val(args)
    path index
    path catalogue
    path exclude
    path species
    path code, stageAs: 'code*/*'

    output:
    tuple val(name), path("species_${name}"), emit: index

    script:
    """
    ${params.kfp} species-index ${index} species_${name} --catalogue ${catalogue} \\
        --exclude ${exclude} --species ${species} ${args}
    """

    stub:
    "mkdir species_${name}"
}

process SPECIES_FIT {
    tag "sample ${sample} ${arm}"
    label 'process_single'
    publishDir "${params.outdir}/species", mode: 'copy', saveAs: { f -> "sample${sample}_${arm}_${f}" }

    input:
    tuple val(sample), path(profile), val(arm), val(args), path(species_index)
    path code, stageAs: 'code*/*'

    output:
    tuple val(sample), val(arm), path('species.tsv'), path('presence.tsv'), path('pfam_presence.tsv'), path('species_units.tsv'), path('function_species.tsv'), emit: fit
    path 'summary.json'

    script:
    // pfam_presence.tsv and function_species.tsv need Pfam labels and hits_em: empty otherwise
    """
    ${params.kfp} species ${profile} ${species_index} species.tsv --units species_units.tsv \\
        --presence presence.tsv --pfam-presence pfam_presence.tsv \\
        --function-taxon function_species.tsv --summary summary.json ${params.fit_args} ${args}
    touch pfam_presence.tsv function_species.tsv
    """

    stub:
    "touch species.tsv presence.tsv pfam_presence.tsv species_units.tsv function_species.tsv summary.json"
}

process REPS {
    label 'process_single'

    input:
    path catalogue
    path species
    path code, stageAs: 'code*/*'

    output:
    path 'reps', emit: reps

    script:
    "${params.bench} uhgg-reps --catalogue ${catalogue} --species ${species} --out-dir reps"

    stub:
    "mkdir reps"
}

process ANNOTATE_REPS {
    label 'process_medium'

    input:
    path index
    path reps
    path catalogue  // the representatives' proteins
    path code, stageAs: 'code*/*'

    output:
    path 'genome_index', emit: index

    script:
    "${params.kfp} annotate-genomes ${index} ${reps}/genomes.tsv genome_index"

    stub:
    "mkdir genome_index"
}

process REPS_SPECIES_INDEX {
    label 'process_medium'

    input:
    path index
    path genome_index
    path code, stageAs: 'code*/*'

    output:
    path 'reps_species', emit: index

    script:
    // kfp-prior's carriage table: one genome per species, shrunk towards genus and family
    "${params.kfp} species-index ${index} reps_species --genomes ${genome_index}"

    stub:
    "mkdir reps_species"
}

process GENOME_FIT {
    tag "sample ${sample}"
    label 'process_single'
    publishDir "${params.outdir}/genomes", mode: 'copy', saveAs: { f -> "sample${sample}_${f}" }

    input:
    tuple val(sample), path(profile)
    path genome_index
    path carriage
    path code, stageAs: 'code*/*'

    output:
    tuple val(sample), path('genomes.tsv'), path('function_taxon.tsv'), emit: fit
    tuple val(sample), path('presence.tsv'), path('pfam_presence.tsv'), emit: prior

    script:
    """
    ${params.kfp} genomes ${profile} ${genome_index} genomes.tsv --function-taxon function_taxon.tsv \\
        ${params.genome_args}
    ${params.kfp_prior} update ${profile} genomes.tsv ${genome_index} ${carriage} \\
        --out presence.tsv --pfam-out pfam_presence.tsv
    touch function_taxon.tsv pfam_presence.tsv
    """

    stub:
    "touch genomes.tsv function_taxon.tsv presence.tsv pfam_presence.tsv"
}

process SYLPH_DB {
    label 'process_medium'
    container params.sylph_container

    input:
    path catalogue

    output:
    path 'reps.syldb', emit: db

    script:
    "sylph sketch -t ${task.cpus} -c ${params.sylph_c} -o reps -g ${catalogue}/species_catalogue/*/*/genome/*.fna"

    stub:
    "touch reps.syldb"
}

process SYLPH {
    tag "sample ${sample}"
    label 'process_medium'
    container params.sylph_container
    publishDir "${params.outdir}/sylph", mode: 'copy', saveAs: { "sample${sample}.tsv" }

    input:
    tuple val(sample), path(r1), path(r2)
    path db

    output:
    tuple val(sample), path('sylph.tsv'), emit: profile

    script:
    """
    sylph sketch -t ${task.cpus} -c ${params.sylph_c} -1 ${r1} -2 ${r2} -d reads
    sylph profile -t ${task.cpus} ${db} reads/*.sylsp > sylph.tsv
    """

    stub:
    "touch sylph.tsv"
}

process SCORE {
    tag "sample ${sample} ${arm}"
    label 'process_single'

    input:
    tuple val(sample), val(arm), val(kind), path(pred, stageAs: 'pred/*'), path(presence, stageAs: 'presence/*'), path(pfam, stageAs: 'pfam/*'), path(units, stageAs: 'units/*'), path(ft, stageAs: 'ft/*')
    path samples
    path species_index
    path code, stageAs: 'code*/*'

    output:
    path 'species_score.tsv', emit: score
    path 'species_calibration.tsv', emit: calibration

    script:
    """
    ${params.bench} uhgg-score --truth ${samples} --species-index ${species_index} --sample ${sample} \\
        --arm ${arm} --kind ${kind} ${opt('pred', pred)} ${opt('presence', presence)} \\
        ${opt('pfam-presence', pfam)} ${opt('units', units)} ${opt('function-taxon', ft)}
    """

    stub:
    "touch species_score.tsv species_calibration.tsv"
}

workflow {
    def bench_py = file("${projectDir}/../fmh-benchmark/bench.py")
    // ponytail: tracks the Python package only; Rust kernel changes still need a fresh run
    def ch_code = channel.fromPath("${projectDir}/../../python/{kmer_functional_profiler,kfp_prior}/*.py").collect()
    METADATA()
    PICK(METADATA.out.metadata, bench_py)
    FETCH_CATALOGUE(METADATA.out.metadata, PICK.out.species)
    FETCH_GENOMES(PICK.out.samples)
    SIMULATE(channel.of(1..params.replicates), PICK.out.samples, FETCH_GENOMES.out.genomes, bench_py)
    def ch_index = params.index
        ? channel.value(file(params.index, checkIfExists: true))
        : INDEX(file(params.members, checkIfExists: true), params.pfam ? file(params.pfam, checkIfExists: true) : [], ch_code).index
    PROFILE(SIMULATE.out.reads, ch_index, ch_code)

    // Species model and its ablations: one species index per variant, arms on top
    SPECIES_INDEX(
        channel.fromList(params.species_indexes).map { v -> [v.name, v.args] },
        ch_index, FETCH_CATALOGUE.out.catalogue, PICK.out.exclude, PICK.out.species, ch_code,
    )
    def ch_arms = channel.fromList(params.species_arms).map { a -> [a.index, a.arm, a.args] }
        .combine(SPECIES_INDEX.out.index, by: 0)  // index, arm, args, species index
    SPECIES_FIT(
        PROFILE.out.profile.combine(ch_arms)
            .map { sample, profile, _index, arm, args, si -> [sample, profile, arm, args, si] },
        ch_code,
    )
    def ch_default = SPECIES_INDEX.out.index.filter { it[0] == 'default' }.map { it[1] }.first()

    // Genome mode on the representatives (G1), kfp-prior (B1) on its fit; sylph
    REPS(FETCH_CATALOGUE.out.catalogue, PICK.out.species, bench_py)
    ANNOTATE_REPS(ch_index, REPS.out.reps, FETCH_CATALOGUE.out.catalogue, ch_code)
    REPS_SPECIES_INDEX(ch_index, ANNOTATE_REPS.out.index, ch_code)
    GENOME_FIT(PROFILE.out.profile, ANNOTATE_REPS.out.index, REPS_SPECIES_INDEX.out.index, ch_code)
    SYLPH(SIMULATE.out.reads, SYLPH_DB(FETCH_CATALOGUE.out.catalogue).db)

    SCORE(
        SPECIES_FIT.out.fit
            .map { sample, arm, sp, pr, pf, units, ft -> [sample, arm, 'species', sp, pr, pf, units, ft] }
            .mix(GENOME_FIT.out.fit.map { sample, g, ft -> [sample, 'genomes', 'genomes', g, [], [], [], ft] })
            .mix(GENOME_FIT.out.prior.map { sample, pr, pf -> [sample, 'kfp_prior', 'species', [], pr, pf, [], []] })
            .mix(SYLPH.out.profile.map { sample, p -> [sample, 'sylph', 'sylph', p, [], [], [], []] }),
        PICK.out.samples,
        ch_default,
        bench_py,
    )
    SCORE.out.score.collectFile(name: 'species_scores.tsv', keepHeader: true, sort: true, storeDir: params.outdir)
    SCORE.out.calibration.collectFile(name: 'species_calibration.tsv', keepHeader: true, sort: true, storeDir: params.outdir)
}
