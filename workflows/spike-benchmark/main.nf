// Spike-in benchmark (plan: phase 13, steps 3-6 and its gate; phase 14, step 5): strains of
// known identity and depth, and withheld siblings, spiked into reads (optionally into a real
// background), then profiled by kfp and followed up by workflows/pileup's PILEUP (DIAMOND
// for amino-acid indexes, bwa-mem2 for nucleotide ones). Nucleotide runs add KMA and kfp's
// amino-acid query on the translated targets. Scored with exact truth by spike.py. See
// README.md.

include { PILEUP } from '../pileup/main.nf'
include { QUERY as QUERY_AA } from '../pileup/main.nf'

process FAMILIES {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path spike  // spike.py, staged so that editing it reruns the task
    path seed_fasta
    val opts  // options as inputs, not params in the script: Nextflow does not hash params

    output:
    tuple path('targets.fa'), path('groups.tsv'), emit: targets

    script:
    """
    ${params.venv}/bin/python ${spike} families --seed-fasta ${seed_fasta} \\
        ${opts}
    """
}

process SIBLINGS {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path spike  // spike.py, staged so that editing it reruns the task
    tuple path(targets), path(groups)
    val opts

    output:
    tuple path('index_targets.fa'), path('index_groups.tsv'), emit: targets
    path 'siblings.tsv', emit: siblings
    path 'withhold.txt', emit: withhold

    script:
    def g = groups.name == 'NO_GROUPS' ? '' : "--groups ${groups}"
    """
    ${params.venv}/bin/python ${spike} siblings --targets ${targets} ${g} \
        ${opts}
    """
}

process SIMULATE {
    tag "seed ${seed}"
    label 'process_low'

    input:
    path spike  // spike.py, staged so that editing it reruns the task
    tuple val(seed), path(targets), path(groups), path(background)
    path siblings
    val opts  // --molecule and simulate_args; --background-pairs with a background

    output:
    tuple val("seed${seed}"), path('reads_1.fq.gz'), path('reads_2.fq.gz'), emit: reads
    tuple val("seed${seed}"), path('truth.tsv'), emit: truth

    script:
    def g = "--groups ${groups} --siblings ${siblings}"
    def bg = background ? "--background ${background.join(' ')}" : ''
    """
    ${params.venv}/bin/python ${spike} simulate --targets ${targets} ${g} \\
        --seed ${seed} ${bg} ${opts}
    """
}

process INDEX_TARGETS {
    tag "${name}"
    label 'process_medium'

    input:
    tuple val(name), path(targets), path(groups), val(args), path(model)  // model: or NO_MODEL

    output:
    tuple val(name), path("index_${name}"), emit: index

    script:
    def g = groups.name == 'NO_GROUPS' ? '' : "--groups ${groups}"
    // the aai survival model (e.g. step34q) on the amino-acid index the pileups check
    def attach = model.name == 'NO_MODEL' ? '' : "\n${params.kfp} aai-model index_${name} ${model}"
    """
    ${params.kfp} index ${targets} index_${name} --targets ${g} ${args}${attach}
    """
}

process TRANSLATE {
    label 'process_low'

    input:
    path spike  // spike.py, staged so that editing it reruns the task
    path targets

    output:
    path 'targets.faa', emit: faa

    script:
    "${params.venv}/bin/python ${spike} translate --targets ${targets} --out targets.faa"
}

process KMA {
    tag "${sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/kma:1.6.18--h118bc1c_0'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/kma.res" }, pattern: 'kma.res'

    input:
    tuple val(sample), path(r1), path(r2)
    path targets
    val opts

    output:
    tuple val(sample), path('kma.res'), emit: res

    script:
    """
    kma index -i ${targets} -o db
    kma -ipe ${r1} ${r2} -o kma -t_db db -t ${task.cpus} ${opts}
    """
}

process SCORE {
    tag "${sample}"
    label 'process_low'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/${it}" }

    input:
    path spike  // spike.py, staged so that editing it reruns the task
    tuple val(sample), path(truth), path(profile), path(pileup), path(ref), path(aa_profile, stageAs: "aa/*"), path(kma)
    path groups

    output:
    path 'scores.tsv', emit: scores
    path 'units.tsv'

    script:
    def aa = aa_profile.toString().endsWith('NO_AA') ? '' : "--aa-profile ${aa_profile}"
    def k = kma.name == 'NO_KMA' ? '' : "--kma ${kma}"
    def g = groups.name == 'NO_GROUPS' ? '' : "--groups ${groups}"
    """
    ${params.venv}/bin/python ${spike} score --sample ${sample} --truth ${truth} \\
        --profile ${profile} --pileup ${pileup} --ref ${ref} ${aa} ${k} ${g}
    """
}

process SUMMARY {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path spike  // spike.py, staged so that editing it reruns the task
    path scores, stageAs: 'scores/*/scores.tsv'

    output:
    path 'summary.tsv'

    script:
    "${params.venv}/bin/python ${spike} summary ${scores} --out summary.tsv"
}

// a placeholder input, never read
def none(name) {
    file("${projectDir}/assets/${name}")
}

workflow {
    if (!params.targets && !params.seed_fasta) {
        error 'Set --targets (a FASTA, with optional --groups) or --seed_fasta (and --families) to make families'
    }
    def spike = file("${projectDir}/spike.py")
    def given = params.targets
        ? channel.value(tuple(file(params.targets, checkIfExists: true),
            params.groups ? file(params.groups, checkIfExists: true) : none('NO_GROUPS')))
        : FAMILIES(spike, file(params.seed_fasta, checkIfExists: true),
            "--families ${params.families} --molecule ${params.molecule} ${params.families_args}").targets
    def sibs = SIBLINGS(spike, given, "--molecule ${params.molecule} ${params.siblings_args}")
    def targets = sibs.targets.first()  // the targets and their siblings, every target's unit
    def background = params.background
        ? channel.value(params.background.toString().tokenize(',').collect { file(it, checkIfExists: true) })
        : channel.value([])
    def sim = SIMULATE(spike, channel.of(1..params.seeds).combine(targets).combine(background.map { [it] }),
        sibs.siblings.first(),
        "--molecule ${params.molecule} ${params.simulate_args}" + (
            params.background ? " --background-pairs ${params.background_pairs}" : ''))
    def nucleotide = params.molecule == 'nucleotide'
    def to_index = targets.map { fa, g ->
        tuple('kfp', fa, g, (nucleotide ? '--alphabet dna ' : '') + params.index_args,
            params.aai_model && !nucleotide ? file(params.aai_model, checkIfExists: true) : none('NO_MODEL'))
    }
    if (nucleotide && params.aa_arm) {
        to_index = to_index.mix(TRANSLATE(spike, targets.map { fa, g -> fa }).faa
            .combine(targets.map { fa, g -> g })
            .map { faa, g -> tuple('aa', faa, g, params.aa_index_args, none('NO_MODEL')) })
    }
    def indexes = INDEX_TARGETS(to_index).index.branch { name, dir ->
        kfp: name == 'kfp'
            return dir
        aa: true
            return dir
    }
    def pileup = PILEUP(sim.reads, indexes.kfp.first(), channel.value([]), sibs.withhold.first())
    def aa = nucleotide && params.aa_arm
        ? QUERY_AA(sim.reads.map { s, r1, r2 -> tuple("${s}~aa", r1, r2) }, indexes.aa.first(),
            indexes.aa.first().map { it.resolve('meta.json') }, params.query_args)
            .profile.map { s, p -> tuple(s - '~aa', p) }
        : sim.reads.map { s, r1, r2 -> tuple(s, none('NO_AA')) }
    def kma = nucleotide && params.kma
        ? KMA(sim.reads, targets.map { fa, g -> fa }, params.kma_args).res
        : sim.reads.map { s, r1, r2 -> tuple(s, none('NO_KMA')) }

    def scores = SCORE(
        spike,
        sim.truth.join(pileup.profile).join(pileup.pileup).join(pileup.ref).join(aa).join(kma),
        targets.map { fa, g -> g }
    ).scores
    SUMMARY(spike, scores.collect())
}
