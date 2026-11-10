// The tree index's route on its own, on the compact index (docs/tree-index-plan.md, Progress log:
// the compact route, 2026-11-10): tier 1 on MGnify90 units, confirmation on each unit's
// representative, every read pair recruited to the confirmed representatives by seed votes and
// pseudo-alignment, scored by depth x AAI. Each process runs one step of
// ../pfam-anchor/route90.py (or ../tree-index/compact.py) on the shared --data directory and
// passes a token. On HPC: hpc/kfp-ablations/route90/run.sh.
//
// HOLDOUT  held-out proteins and their nearest relative left in the reference   recall.py holdout
// READS    strains and relatives at depths 0.5-30x                              route90.py reads
// SHARD    one shard's FASTA, unit offsets and keys (k 9, sqrt(members)), per    compact.py shard
//          index size: 1 in 32 (standard), 1 in 8 and 1 in 4 (the larger options)
// BUILD    the shards' keys merged, per index size (--rates)                    route90.py build
// QUERY    tier 1 and confirm (compact.py), then recruit (recruit.py, apart)     route90.py query
// SCORE    route90_recall_r<rate>.tsv, route90_targets_r<rate>.tsv              route90.py score

nextflow.enable.dsl = 2

def python() { params.python ?: "${params.venv}/bin/python" }
def data() { params.data ?: "${file(params.outdir).toAbsolutePath()}/data" }

// route90.py <step> on the shared data directory for one index size (--rate), threads and the
// watchdog's cap from the task
def step(name, task, rate) {
    """
    export RAYON_NUM_THREADS=${task.cpus} POLARS_MAX_THREADS=${task.cpus}
    ${python()} ${projectDir}/../pfam-anchor/route90.py ${name} --out ${data()} \
        --holdout ${data()}/holdout --max-gb ${(task.memory.toGiga() * 0.9).round(1)} \
        --strains ${params.strains} --shards ${params.shards} --rate ${rate}
    """
}

process HOLDOUT {
    label 'process_medium'
    output: val 'holdout'
    script:
    "${python()} ${projectDir}/../tree-index/recall.py holdout --members ${params.members} --per-bin ${params.per_bin} --pool ${params.pool} --seed 1 --out ${data()}/holdout"
}

process READS {
    label 'process_medium'
    input: val token
    output: val 'reads'
    script: step('reads', task, 32)
}

process SHARD {
    label 'process_medium'
    input: tuple val(token), val(rate), val(shard)
    output: val rate
    script:
    """
    export RAYON_NUM_THREADS=${task.cpus} POLARS_MAX_THREADS=${task.cpus}
    ${python()} ${projectDir}/../tree-index/compact.py shard --members ${data()}/holdout/members.parquet \
        --index ${data()}/index-r${rate} --shards ${params.shards} --shard ${shard} --rate ${rate}
    """
}

process BUILD {
    label 'process_medium'
    input: val rate
    output: val rate
    script: step('build', task, rate)
}

process QUERY {
    label 'process_high'
    input: val rate
    output: val rate
    script: step('query', task, rate)
}

process SCORE {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'
    input: val rate
    output: path '*.{tsv,json}'
    script:
    step('score', task, rate) + "cp ${data()}/route90_*_r${rate}.tsv ${data()}/*_times_r${rate}.json .\n"
}

workflow {
    def n = params.shards as int
    held = HOLDOUT()
    rates = Channel.fromList(params.rates.toString().tokenize(',')*.trim())
    built = SHARD(held.combine(rates).combine(Channel.of(0..<n)))
        .map { [it, 1] }.groupTuple(size: n).map { it[0] }  // a rate once all its shards are built
    SCORE(QUERY(BUILD(READS(held).combine(built).map { it[1] })))
}
