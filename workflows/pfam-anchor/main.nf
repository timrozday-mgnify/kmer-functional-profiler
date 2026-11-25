// The Pfam-anchored route and its benchmark (docs/pfam-anchor-plan.md), end to end on one
// reference. On HPC: hpc/kfp-ablations/pfam-anchor/run.sh. Locally: -profile test,laptop.
//
// FIXTURE     test profile only: pfam-reannotate's fixture -> members.parquet, families.tsv
// ANNOTATE    the families' HMMs on every member (hmmsearch)              prepare.py annotate
// NULLREF     the members reversed, the bar on weak reference domains      prepare.py nullref
// READS       held-out relatives, simulated strains and their reads        prepare.py reads
// SEARCH      six-frame ORFs and their hits, and the reversed reads'       prepare.py search
// RECRUIT     kfpfam's weak-read recruitment (--recruit)                   prepare.py recruit
// DIAMOND     every mate against the reference, the method-free yardstick  prepare.py diamond
// REPORT      hit recall by identity bin                                   prepare.py report
// ORACLE      every read hit against every domain of its family            anchor.py oracle
// ANCHORED    the anchored route (anchor residues, column identity)        anchor.py query
// COMPARE     the route against DIAMOND                                    anchor.py compare
// UNANCHORED  tree index tiers 1 and 2 on the whole reference              anchor.py unanchored
// CONFIRM     anchored and unanchored units -> the confirmed reference     anchor.py confirm
// TREE        the tree index on the confirmed reference; tier 2 alone      anchor.py treeindex
// SEEDLINES   tier 2's seed alignments as call's lines; kfp-pileup call    anchor.py seedlines
// MAPALL      DIAMOND on the confirmed reference; kfp-pileup call          anchor.py mapall
// SCORE       mapping, gather and pileups against the truth                anchor.py mapcmp, truthscore
// BENCH       phase 6: recall by depth x AAI, misses by stage, cost        bench.py score, cost
//
// Every step reads and writes --data (the scripts' fixed layout), so processes pass a token,
// not files; outside tools run under --runtime (docker, singularity, or none).

nextflow.enable.dsl = 2

params.members  = ''        // reference members parquet: protein_id, cluster_rep, full_length, sequence
params.families = ''        // pfam_accession (and proteins) TSV: the families searched
params.pfam     = ''        // Pfam-A.hmm[.gz], the HMMs to cut the families from
params.kfpfam   = ''        // a kfpfam checkout (recruit/, data/prep/), for --recruit
params.recruit  = true      // kfpfam's recruited read hits; false: every hmmsearch hit
params.strains  = 2000      // simulated strains (prepare.py reads)
params.per_bin  = 250       // held-out relatives per identity bin
params.pool     = 15000     // proteins searched for relatives
params.min_weight = 0.1     // confirm: called pairs per unit
params.depths   = '0.5 1 2 5 10 30'  // phase 6's depth grid; '' = the earlier lognormal
params.background_1 = ''    // real read pairs (R1, R2) added to the strains set
params.background_2 = ''
params.background_pairs = 2000000
params.annotated = ''       // another run's data directory: reuse its annotate and nullref
params.runtime  = 'singularity'
params.outdir   = 'results'
params.data     = ''        // the scripts' shared data directory; '' = <outdir>/data
params.venv     = "${projectDir}/../../.venv"
params.release  = ''        // test profile: pfam-reannotate's fixture tables

def py() { "${params.venv}/bin/python" }
def data() { file(params.data ?: "${params.outdir}/data").toAbsolutePath() }
def prep(step, task, extra = '') {
    """${py()} ${projectDir}/prepare.py ${step} --out ${data()} --members ${file(params.members).toAbsolutePath()} \\
        --families ${file(params.families).toAbsolutePath()} --pfam ${file(params.pfam).toAbsolutePath()} \\
        --strains ${params.strains} --per-bin ${params.per_bin} --pool ${params.pool} \\
        --threads ${task.cpus} --max-gb ${task.memory.toGiga() * 0.9} --runtime ${params.runtime} \\
        ${params.depths ? "--depths ${params.depths}" : ''} \\
        ${params.background_1 ? "--background ${file(params.background_1).toAbsolutePath()} ${file(params.background_2).toAbsolutePath()} --background-pairs ${params.background_pairs}" : ''} ${extra}"""
}
def hits() { params.recruit ? 'recruited_top' : 'hmmsearch' }
def anch(step, task, extra = '') {
    """RAYON_NUM_THREADS=${task.cpus} POLARS_MAX_THREADS=${task.cpus} ${py()} ${projectDir}/anchor.py ${step} \\
        --data ${data()} --threads ${task.cpus} --max-gb ${task.memory.toGiga() * 0.9} --runtime ${params.runtime} \\
        --hits ${hits()} --ref weak --min-weight ${params.min_weight} ${extra}"""
}

process FIXTURE {
    label 'process_low'
    input: val ready
    output: val ready
    script:
    """
    ${py()} ${projectDir}/prepare.py fixture --release ${file(params.release).toAbsolutePath()} \\
        --pfam ${file(params.pfam).toAbsolutePath()} --members ${file(params.members).toAbsolutePath()} \\
        --families ${file(params.families).toAbsolutePath()}
    """
}

def reuse(files) {
    def src = file(params.annotated).toAbsolutePath()
    "mkdir -p ${data()} && " + files.collect { "cp ${src}/${it} ${data()}/" }.join(' && ')
}

process ANNOTATE {
    label 'process_medium'
    input: val ready
    output: val ready
    script: params.annotated ? reuse(['hmms.hmm', 'ga.tsv', 'hits.parquet', 'hits.stats.tsv']) : prep('annotate', task)
}

process NULLREF {
    label 'process_medium'
    input: val ready
    output: val ready
    script: params.annotated ? reuse(['null_hits.parquet', 'null_hits.stats.tsv']) : prep('nullref', task)
}

process READS {
    label 'process_medium'
    input: val ready
    output: val ready
    script: prep('reads', task)
}

process SEARCH {
    label 'process_medium'
    input: val ready
    output: val ready
    script: prep('search', task)
}

process RECRUIT {
    label 'process_medium'
    input: val ready
    output: val ready
    script: prep('recruit', task, "--kfpfam ${file(params.kfpfam).toAbsolutePath()}")
}

process DIAMOND {
    label 'process_medium'
    input: val ready
    output: val ready
    script: prep('diamond', task)
}

process REPORT {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'
    input: val ready
    output: path 'report.tsv'
    script: "${prep('report', task)} && cp ${data()}/report.tsv ."
}

process ORACLE {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'
    input: val ready
    output: tuple val(ready), path('oracle.tsv'), path('alignment.tsv')
    script: "${anch('oracle', task)} && cp ${data()}/oracle.tsv ${data()}/alignment.tsv ."
}

process ANCHORED {
    label 'process_medium'
    input: val ready
    output: val ready
    script: anch('query', task)
}

process COMPARE {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'
    input: val ready
    output: path 'phase3*.tsv'
    script: "${anch('compare', task)} && cp ${data()}/phase3_${hits()}_weak.tsv ${data()}/status_${hits()}_weak.tsv ."
}

process UNANCHORED {
    label 'process_high'
    input: val ready
    output: val ready
    script: anch('unanchored', task)
}

process CONFIRM {
    label 'process_low'
    input: val ready
    output: val ready
    script: anch('confirm', task, '--with-unanchored')
}

process TREE {
    label 'process_medium'
    input: val ready
    output: val ready
    script: anch('treeindex', task, '--with-unanchored')
}

process SEEDLINES {
    label 'process_medium'
    input: val ready
    output: val ready
    script: anch('seedlines', task, '--with-unanchored')
}

process MAPALL {
    label 'process_medium'
    input: val ready
    output: val ready
    script: anch('mapall', task, '--with-unanchored')
}

process SCORE {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'
    input: val ready
    output: path '*.tsv'
    script:
    def tag = "confirmed_w${params.min_weight}_u".replace('.0_u', '_u')
    """
    ${anch('mapcmp', task, '--with-unanchored')}
    ${anch('truthscore', task, '--with-unanchored')}
    cp ${data()}/mapall_${tag}.tsv ${data()}/truthscore_${tag}.tsv .
    """
}

process BENCH {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'
    input: val ready
    output: path 'bench_*'
    script:
    def args = "--data ${data()} --hits ${hits()} --ref weak --min-weight ${params.min_weight}"
    """
    ${py()} ${projectDir}/bench.py score ${args}
    ${py()} ${projectDir}/bench.py cost ${args} --trace ${file(params.outdir).toAbsolutePath()}/trace.tsv
    cp ${data()}/bench_*.tsv .
    """
}

workflow {
    // every process passes this token on: the task hash doesn't see params read inside prep()
    // and anch(), so a changed parameter must reach every task through its input
    token = Channel.value(params.findAll { it.key != 'outdir' }.toString().md5())
    start = params.release ? FIXTURE(token) : token
    annotated = ANNOTATE(start)
    nullref = NULLREF(annotated)
    reads = READS(annotated)
    searched = SEARCH(reads)
    hit_ready = params.recruit ? RECRUIT(searched) : searched
    dia = DIAMOND(searched)
    REPORT(dia)
    oracle = ORACLE(dia.combine(hit_ready).map { it[0] })
    anchored = ANCHORED(oracle.map { it[0] }.combine(nullref).map { it[0] })
    COMPARE(anchored)
    unanchored = UNANCHORED(reads)
    confirmed = CONFIRM(anchored.combine(unanchored).map { it[0] })
    tree = TREE(confirmed)
    seeded = SEEDLINES(tree)
    mapped = MAPALL(confirmed)
    scored = SCORE(seeded.combine(mapped).map { it[0] })
    BENCH(scored.map { it[0] })
}
