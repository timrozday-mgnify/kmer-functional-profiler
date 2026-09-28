// MGnify Proteins subset (whole MGnify90 clusters of one biome) -> members table -> index.
// See README.md.

process MEMBERSHIP {
    label 'process_medium'

    output:
    path 'membership.parquet', emit: membership

    script:
    """
    ${params.python} ${projectDir}/mgnify_subset.py --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' \\
        membership --release '${params.release}' --biome '${params.biome}' --sample ${params.sample}
    """

    stub:
    "touch membership.parquet"
}

process EXTRACT {
    tag "shard ${shard}"
    label 'process_medium'

    input:
    tuple val(shard), val(lo), val(hi)
    path membership

    output:
    path "shard${shard}.*.parquet", emit: parts

    script:
    """
    ${params.python} ${projectDir}/mgnify_subset.py --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' \\
        extract --release '${params.release}' --membership ${membership} \\
        --lo ${lo} --hi ${hi} --prefix shard${shard}
    """

    stub:
    "touch shard${shard}.seqs.parquet shard${shard}.pfam.parquet"
}

process MERGE {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'

    input:
    path membership
    path parts

    output:
    path 'members.parquet', emit: members
    path 'pfam.parquet', emit: pfam

    script:
    def prefixes = parts.collect { it.name.tokenize('.')[0] }.unique().sort().join(' ')
    """
    ${params.python} ${projectDir}/mgnify_subset.py --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' \\
        merge --membership ${membership} ${prefixes}
    """

    stub:
    "touch members.parquet pfam.parquet"
}

process INDEX {
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy'

    input:
    path members
    path pfam

    output:
    path 'index', emit: index

    script:
    """
    ${params.kfp} index ${members} index --pfam ${pfam} ${params.index_args}
    """

    stub:
    "mkdir index && touch index/meta.json"
}

workflow {
    MEMBERSHIP()
    // Even protein_id ranges; MERGE fails if max_protein_id leaves members out.
    def width = (params.max_protein_id as long).intdiv(params.shards as int) + 1
    ch_shards = channel.of(0..<(params.shards as int)).map { i -> [i, i * width, (i + 1) * width] }
    EXTRACT(ch_shards, MEMBERSHIP.out.membership)
    MERGE(MEMBERSHIP.out.membership, EXTRACT.out.parts.flatten().collect())
    if (params.index) {
        INDEX(MERGE.out.members, MERGE.out.pfam)
    }
}
