// Presence calibration at full-index scale (plan, phase 13 / phase 7: presence at a false
// discovery rate): PICK the index's own units, SIMULATE strains of them into a real background
// (spike-benchmark's spike.py), REVERSE the background into a null sample, QUERY every sample
// at each --background-share, SCORE the spiked samples against the unspiked and null ones,
// DECIDE. See README.md; on HPC, hpc/kfp-ablations/presence-calibration/run.sh.

include { INDEX } from '../pileup/main.nf'

def python() {
    "${params.venv}/bin/python"
}

process PICK {
    tag "seed ${seed}"
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', pattern: 'units.tsv', saveAs: { "seed${seed}/${it}" }

    input:
    path code
    val seed
    path index
    path meta, stageAs: 'index.meta.json'  // a rebuilt index reruns the pick
    val store

    output:
    tuple val(seed), path('targets.fa'), path('groups.tsv'), path('units.tsv'), emit: picked

    script:
    def s = store ? "--store ${store}" : ''
    """
    ${python()} ${code} pick --index ${index} ${s} --units ${params.units} \\
        --members ${params.members} --seed ${seed} --out .
    """
}

process SIMULATE {
    tag "seed ${seed}"
    label 'process_low'
    publishDir params.outdir, mode: 'copy', pattern: 'truth.tsv', saveAs: { "seed${seed}/${it}" }

    input:
    path spike
    tuple val(seed), path(targets), path(groups), path(units), path(bg1), path(bg2)

    output:
    tuple val("seed${seed}"), path('reads_1.fq.gz'), path('reads_2.fq.gz'), emit: reads
    tuple val("seed${seed}"), path('truth.tsv'), path(units), emit: truth

    script:
    """
    ${python()} ${spike} simulate --targets ${targets} --groups ${groups} --seed ${seed} \\
        --present 1 ${params.simulate_args} --background ${bg1} ${bg2} --out .
    """
}

process REVERSE {
    label 'process_low'

    input:
    path code
    tuple path(r1), path(r2)

    output:
    tuple val('null'), path('null_1.fq.gz'), path('null_2.fq.gz'), emit: reads

    script:
    "${python()} ${code} reverse --reads ${r1} ${r2} --out null_1.fq.gz null_2.fq.gz"
}

process QUERY {
    tag "${sample} share ${share}"
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/share${share}/${it}" }

    input:
    tuple val(sample), path(r1), path(r2), val(share)
    path index
    path meta, stageAs: 'index.meta.json'
    path code, stageAs: 'code/*'  // unused: staged so that -resume reruns on query changes

    output:
    tuple val(sample), val(share), path('profile.tsv'), emit: profile
    path 'stats.json'

    script:
    """
    ${params.kfp} query ${index} ${r1} ${r2} --out profile.tsv --stats stats.json \\
        --background-share ${share} ${"${params.aai}" == 'true' ? '--aai' : ''} \\
        --gather ${params.gather} ${params.query_args}
    """
}

process SCORE {
    tag "${sample} share ${share}"
    label 'process_low'

    input:
    path code
    tuple val(sample), val(share), path(profile), path(truth), path(units), path(unspiked, stageAs: 'unspiked/*'), path(nul, stageAs: 'null/*')

    output:
    path 'scores.tsv', emit: scores

    script:
    """
    ${python()} ${code} score --sample ${sample} --share ${share} --truth ${truth} \\
        --units ${units} --profile ${profile} --unspiked ${unspiked} --null ${nul} --out scores.tsv
    """
}

process DECIDE {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path code
    path scores, stageAs: 'scores/*.tsv'

    output:
    path 'summary.tsv'
    path 'decision.tsv'
    path 'scores.tsv'

    script:
    """
    ${python()} ${code} decide ${scores} --q ${params.q} --max-false ${params.max_false} \\
        --max-null ${params.max_null}
    """
}

// one sample's profiles keyed by share
def by_share(profiles, name) {
    profiles.filter { it[0] == name }.map { s, share, p -> tuple(share, p) }
}

workflow {
    if (!params.background) {
        error 'Set --background r1.fq.gz,r2.fq.gz (the real reads the spikes go into)'
    }
    if (!params.index && !params.targets) {
        error 'Set --index (a kfp index directory, absolute) or --targets (a FASTA to index)'
    }
    def code = file("${projectDir}/calibrate.py", checkIfExists: true)
    def spike = file("${projectDir}/../spike-benchmark/spike.py", checkIfExists: true)
    def query_code = ['query', 'partition'].collect {
        file("${projectDir}/../../python/kmer_functional_profiler/${it}.py", checkIfExists: true)
    }
    def index = params.index
        ? channel.value(file(params.index, checkIfExists: true))
        : INDEX(file(params.targets, checkIfExists: true), params.index_args).index
    def meta = index.map { it.resolve('meta.json') }
    def bg = params.background.tokenize(',').collect { file(it, checkIfExists: true) }
    def bg_pair = channel.value(tuple(bg[0], bg[1]))
    def store = params.store ? files(params.store, checkIfExists: true).collect { it.toString() }.join(' ') : ''

    def picked = PICK(code, channel.of(1..params.seeds), index, meta, store).picked
    def sim = SIMULATE(spike, picked.combine(bg_pair))
    def unspiked = channel.value(tuple('unspiked', bg[0], bg[1]))
    def samples = sim.reads.mix(unspiked, REVERSE(code, bg_pair).reads)
    def shares = channel.fromList(params.shares.toString().tokenize(','))
    def profiles = QUERY(samples.combine(shares), index, meta, query_code).profile

    // each spiked sample's profile beside the unspiked and null profiles at its share
    def spiked = profiles.filter { it[0].startsWith('seed') }
        .combine(sim.truth, by: 0)  // sample, share, profile, truth, units
        .map { s, share, p, t, u -> tuple(share, s, p, t, u) }
        .combine(by_share(profiles, 'unspiked'), by: 0)
        .combine(by_share(profiles, 'null'), by: 0)
        .map { share, s, p, t, u, up, np -> tuple(s, share, p, t, u, up, np) }
    DECIDE(code, SCORE(code, spiked).scores.collect())
}
