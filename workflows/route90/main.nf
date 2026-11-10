// The tree index's route on its own, on the compact index (docs/tree-index-plan.md, Progress log:
// the compact route, 2026-11-10): tier 1 on MGnify90 units, confirmation on each unit's
// representative, every read pair recruited to the confirmed representatives by seed votes and
// pseudo-alignment, scored by depth x AAI. Each process runs one step of
// ../pfam-anchor/route90.py (or ../tree-index/compact.py) on the shared --data directory and
// passes a token. On HPC: hpc/kfp-ablations/route90/run.sh.
//
// HOLDOUT  held-out proteins and their nearest relative left in the reference   recall.py holdout
// READS    strains and relatives at depths 0.5-30x                              route90.py reads
// SHARD    one shard's FASTA, unit offsets and keys (k 9, 1 in 32, sqrt)        compact.py shard
// BUILD    the shards' keys merged                                              route90.py build
// QUERY    tier 1 and confirm (compact.py), then recruit (recruit.py, apart)     route90.py query
// SCORE    route90_recall.tsv, route90_targets.tsv                              route90.py score

nextflow.enable.dsl = 2

def python() { params.python ?: "${params.venv}/bin/python" }
def data() { params.data ?: "${file(params.outdir).toAbsolutePath()}/data" }

// route90.py <step> on the shared data directory, threads and the watchdog's cap from the task
def step(name, task) {
    """
    export RAYON_NUM_THREADS=${task.cpus} POLARS_MAX_THREADS=${task.cpus}
    ${python()} ${projectDir}/../pfam-anchor/route90.py ${name} --out ${data()} \
        --holdout ${data()}/holdout --max-gb ${(task.memory.toGiga() * 0.9).round(1)} \
        --strains ${params.strains} --shards ${params.shards}
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
    script: step('reads', task)
}

process SHARD {
    label 'process_medium'
    input: tuple val(token), val(shard)
    output: val 'shard'
    script:
    """
    export RAYON_NUM_THREADS=${task.cpus} POLARS_MAX_THREADS=${task.cpus}
    ${python()} ${projectDir}/../tree-index/compact.py shard --members ${data()}/holdout/members.parquet \
        --index ${data()}/index --shards ${params.shards} --shard ${shard}
    """
}

process BUILD {
    label 'process_medium'
    input: val token
    output: val 'build'
    script: step('build', task)
}

process QUERY {
    label 'process_high'
    input: val token
    output: val 'query'
    script: step('query', task)
}

process SCORE {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'
    input: val token
    output: path '*.{tsv,json}'
    script:
    step('score', task) + "cp ${data()}/route90_*.tsv ${data()}/*_times.json .\n"
}

workflow {
    held = HOLDOUT()
    shards = SHARD(held.combine(Channel.of(0..<(params.shards as int))))
    SCORE(QUERY(BUILD(READS(held).combine(shards.collect()).map { it[0] })))
}
