// MGnify Proteins subset (whole MGnify90 clusters of one biome) -> members table -> index,
// the full-scale cost study (per-cluster statistics -> predicted index size), and the
// partitioned index build over the buckets.
// See README.md.

process MEMBERSHIP {
    label 'process_medium'

    input:
    val sample
    val width

    output:
    path 'membership', emit: membership

    script:
    """
    ${params.python} ${projectDir}/mgnify_subset.py --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' \\
        membership --release '${params.release}' --biome '${params.biome}' --sample ${sample} \\
        --max-protein-id ${params.max_protein_id} --shard-width ${width}
    """

    stub:
    "mkdir membership"
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
        --lo ${lo} --hi ${hi} --prefix shard${shard} --buckets ${params.buckets} \\
        ${params.pfam ? '--pfam' : '--no-pfam'}
    """

    stub:
    "touch shard${shard}.seqs.parquet shard${shard}.pfam.parquet"
}

process MERGE {
    tag "bucket ${bucket}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_mode

    input:
    val bucket
    path parts

    output:
    path 'members*.parquet', emit: members
    path 'pfam.parquet', emit: pfam, optional: true

    script:
    def prefixes = parts.collect { it.name.tokenize('.')[0] }.unique().sort().join(' ')
    def members = (params.buckets as int) > 1 ? "members.bucket${bucket}.parquet" : 'members.parquet'
    def pfam = params.pfam && bucket == 0 ? '--pfam pfam.parquet' : ''
    """
    ${params.python} ${projectDir}/mgnify_subset.py --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' \\
        merge --bucket ${bucket} --members ${members} ${pfam} ${prefixes}
    """

    stub:
    "touch members.parquet pfam.parquet"
}

process SUBSET {
    tag "1 in ${sample}"
    label 'process_medium'

    input:
    val sample
    path members
    path pfam

    output:
    tuple val(sample), path("members.1in${sample}.parquet"), path("pfam.1in${sample}.parquet")

    script:
    """
    ${params.python} ${projectDir}/mgnify_subset.py --threads ${task.cpus} --memory '${task.memory.toGiga()}GB' \\
        subset --members ${members} --pfam ${pfam} --sample ${sample} \\
        --out-members members.1in${sample}.parquet --out-pfam pfam.1in${sample}.parquet
    """

    stub:
    "touch members.1in${sample}.parquet pfam.1in${sample}.parquet"
}

process INDEX {
    tag "1 in ${sample}"
    label 'process_high_memory'
    publishDir path: { "${params.outdir}/1in${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(members), path(pfam)

    output:
    path 'index', emit: index

    script:
    """
    ${params.kfp} index ${members} index --pfam ${pfam} ${params.index_args}
    """

    stub:
    "mkdir index && touch index/meta.json"
}

process STATS {
    tag "${members.baseName}"
    label 'process_medium'

    input:
    path members

    output:
    path "${members.baseName}.*.parquet", emit: stats

    script:
    """
    ${params.python} ${projectDir}/mgnify_subset.py stats --members ${members} \\
        --prefix ${members.baseName} --rate ${params.stats_rate} ${params.stats_args}
    """

    stub:
    "touch ${members.baseName}.clusters.parquet ${members.baseName}.pairs.parquet"
}

process COMBINE {
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy'

    input:
    path stats

    output:
    path 'clusters.parquet'
    path 'groups.tsv'
    path 'cost.tsv', emit: cost

    script:
    def prefixes = stats.collect { it.name.replace('.clusters.parquet', '').replace('.pairs.parquet', '') }
        .unique().sort().join(' ')
    def samples = params.cost_samples.toString().tokenize(',').join(' ')
    """
    ${params.python} ${projectDir}/mgnify_subset.py combine ${prefixes} --sample ${samples} \\
        ${params.stats_args} ${params.cost_args}
    """

    stub:
    "touch clusters.parquet groups.tsv cost.tsv"
}

// Partitioned build (kmer_functional_profiler.partition): per bucket, per hash range, once.
process CANDIDATES {
    tag "bucket ${bucket}"
    label 'process_medium'

    input:
    tuple val(bucket), path(members)

    output:
    tuple val(bucket), path(members), path("${members.baseName}.units.parquet"), path("${members.baseName}.stats.json"), emit: units
    path "${members.baseName}.candidates.npy", emit: candidates

    script:
    """
    ${params.python} ${projectDir}/mgnify_subset.py candidates --members ${members} --prefix ${members.baseName} ${params.index_args}
    """

    stub:
    "touch ${members.baseName}.units.parquet ${members.baseName}.stats.json ${members.baseName}.candidates.npy"
}

process BLOOM {
    label 'process_medium'

    input:
    path units
    path candidates

    output:
    path 'bloom.*'

    script:
    def prefixes = candidates.collect { f -> f.name.replace('.candidates.npy', '') }.sort().join(' ')
    """
    ${params.python} ${projectDir}/mgnify_subset.py bloom ${prefixes} --bits-per-key ${params.bloom_bits}
    """

    stub:
    "touch bloom.npy bloom.json"
}

process PRESENCE {
    tag "bucket ${bucket}"
    label 'process_medium'

    input:
    tuple val(bucket), path(members), path(units), path(stats)
    path bloom

    output:
    path "${members.baseName}.range*.parquet"

    script:
    """
    ${params.python} ${projectDir}/mgnify_subset.py presence --members ${members} --prefix ${members.baseName} --bucket ${bucket} \
        --ranges ${params.ranges} ${params.index_args}
    """

    stub:
    "for r in \$(seq 0 ${params.ranges - 1}); do touch ${members.baseName}.range\$r.parquet; done"
}

process GROUPS {
    tag "range ${range}"
    label 'process_medium'

    input:
    tuple val(range), path(presence)

    output:
    path "range${range}.bucket*.parquet", optional: true

    script:
    """
    ${params.python} ${projectDir}/mgnify_subset.py groups --prefix range${range} ${presence}
    """

    stub:
    "touch range${range}.bucket0.parquet"
}

process POSTINGS {
    tag "bucket ${bucket}"
    label 'process_medium'

    input:
    tuple val(bucket), path(members), path(units), path(stats), path(groups)
    path pfam

    output:
    path "${members.baseName}.{final.parquet,final.json,postings.parquet,pfam.parquet}"

    script:
    """
    ${params.python} ${projectDir}/mgnify_subset.py postings --members ${members} --prefix ${members.baseName} --pfam ${pfam} \
        ${params.index_args} ${groups}
    """

    stub:
    "touch ${members.baseName}.final.parquet ${members.baseName}.final.json ${members.baseName}.postings.parquet ${members.baseName}.pfam.parquet"
}

process UNITS {
    label 'process_high_memory'

    input:
    path parts
    path bloom

    output:
    path 'units'

    script:
    def prefixes = parts.findAll { f -> f.name.endsWith('.final.json') }
        .collect { f -> f.name.replace('.final.json', '') }.sort().join(' ')
    """
    ${params.python} ${projectDir}/mgnify_subset.py units ${prefixes} --ranges ${params.ranges} \\
        --pfam ${params.index_args}
    """

    stub:
    "mkdir units && touch units/units.parquet units/pack.json"
}

process PACK_RANGE {
    tag "range ${range}"
    label 'process_medium'

    input:
    tuple val(range), path(postings), path(units)

    output:
    tuple val(range), path("part${range}.*")

    script:
    def prefixes = postings.collect { f -> f.name.replace('.postings.parquet', '') }.sort().join(' ')
    """
    ${params.python} ${projectDir}/mgnify_subset.py pack-range ${prefixes} --range ${range} \\
        --out part${range} ${params.index_args}
    """

    stub:
    "touch part${range}.keys.npy part${range}.json"
}

process CONCAT {
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy'

    input:
    val ranges
    path parts
    path units

    output:
    path 'index'

    script:
    def prefixes = ranges.collect { r -> "part${r}" }.join(' ')
    """
    ${params.python} ${projectDir}/mgnify_subset.py concat ${prefixes} ${params.index_args}
    """

    stub:
    "mkdir index && touch index/meta.json"
}

workflow {
    // Nested subsets: extract at the densest sample, then keep 1 in N of its clusters.
    def samples = params.sample.toString().tokenize(',').collect { it.trim() as long }.sort()
    if (samples.any { it % samples[0] != 0 }) {
        error "--sample ${params.sample}: every sample must be a multiple of the smallest"
    }
    if (params.index && ((params.buckets as int) > 1 || !params.pfam)) {
        error "--index needs --buckets 1 and --pfam true"
    }
    if (params.build && (samples.size() > 1 || !params.pfam)) {
        error "--build needs one --sample and --pfam true"
    }
    // Even protein_id ranges; MEMBERSHIP fails if max_protein_id leaves members out.
    def width = (params.max_protein_id as long).intdiv(params.shards as int) + 1
    MEMBERSHIP(samples[0], width)
    ch_shards = channel.of(0..<(params.shards as int)).map { i -> [i, i * width, (i + 1) * width] }
    EXTRACT(ch_shards, MEMBERSHIP.out.membership)
    MERGE(channel.of(0..<(params.buckets as int)), EXTRACT.out.parts.flatten().collect())
    if (params.stats) {
        STATS(MERGE.out.members)
        COMBINE(STATS.out.stats.flatten().collect())
    }
    if (params.index) {
        SUBSET(channel.fromList(samples), MERGE.out.members.first(), MERGE.out.pfam.first())
        INDEX(SUBSET.out)
    }
    if (params.build) {
        ch_members = MERGE.out.members.flatten().map { f ->
            def m = f.name =~ /bucket(\d+)/
            [m ? m[0][1] as int : 0, f]
        }
        CANDIDATES(ch_members)
        BLOOM(CANDIDATES.out.units.map { t -> t[2] }.collect(), CANDIDATES.out.candidates.collect())
        PRESENCE(CANDIDATES.out.units, BLOOM.out)
        ch_ranges = PRESENCE.out.flatten().map { f -> [(f.name =~ /range(\d+)/)[0][1] as int, f] }
            .groupTuple(size: params.buckets as int)
        GROUPS(ch_ranges)
        ch_groups = GROUPS.out.flatten().map { f -> [(f.name =~ /bucket(\d+)/)[0][1] as int, f] }
            .groupTuple()
        // Buckets with no candidates get no group files.
        ch_postings = CANDIDATES.out.units.join(ch_groups, remainder: true)
            .map { b, members, units, stats, groups -> [b, members, units, stats, groups ?: []] }
        POSTINGS(ch_postings, MERGE.out.pfam.first())
        ch_posted = POSTINGS.out.flatten()
        UNITS(ch_posted.collect(), BLOOM.out)
        ch_range = channel.of(0..<(params.ranges as int))
            .combine(ch_posted.filter { f -> f.name.endsWith('.postings.parquet') }.collect().toList())
            .combine(UNITS.out)
        PACK_RANGE(ch_range)
        ch_parts = PACK_RANGE.out.toSortedList { a, b -> a[0] <=> b[0] }
        CONCAT(
            ch_parts.map { parts -> parts.collect { p -> p[0] } },
            ch_parts.map { parts -> parts.collect { p -> p[1] }.flatten() },
            UNITS.out,
        )
    }
}
