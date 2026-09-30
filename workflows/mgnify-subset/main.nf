// MGnify Proteins subset (whole MGnify90 clusters of one biome) -> members table -> index,
// the full-scale cost study (per-cluster statistics -> predicted index size), the
// partitioned index build over the buckets, and the query cost study (--query).
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

process DEDUP {
    tag "sets ${range}"
    label 'process_medium'

    input:
    val range
    val ranges
    path parts

    output:
    path "sets${range}.*"

    script:
    def prefixes = ranges.collect { r -> "part${r}" }.join(' ')
    """
    ${params.python} ${projectDir}/mgnify_subset.py dedup ${prefixes} --range ${range} \\
        --ranges ${params.ranges} --out sets${range} ${params.index_args}
    """

    stub:
    "touch sets${range}.offsets.npy"
}

process CONCAT {
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy'

    input:
    val ranges
    path parts
    path sets
    path units

    output:
    path 'index'

    script:
    def prefixes = ranges.collect { r -> "part${r}" }.join(' ')
    def sets_ = (0..<(params.ranges as int)).collect { s -> "sets${s}" }.join(' ')
    """
    ${params.python} ${projectDir}/mgnify_subset.py concat ${prefixes} --sets ${sets_} \\
        ${params.index_args}
    """

    stub:
    "mkdir index && touch index/meta.json"
}

// Query cost (phase 6, Q1): nested read subsets of one run x indexes.
process FETCH_READS {
    tag "${run}"

    input:
    val run

    output:
    tuple path("${run}_1.fastq.gz"), path("${run}_2.fastq.gz")

    script:
    """
    read -r _ urls md5s < <(curl -fsS --retry 5 \\
        'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv' \\
        | tail -n 1)
    for url in \${urls//;/ }; do curl -fsSL --retry 5 -C - -O "https://\$url"; done
    paste -d ' ' <(tr ';' '\\n' <<< "\$md5s") <(tr ';' '\\n' <<< "\$urls" | xargs -n 1 basename | sed 's/^/ /') \\
        | md5sum -c -
    """

    stub:
    "touch ${run}_1.fastq.gz ${run}_2.fastq.gz"
}

process LADDER {
    input:
    tuple path(r1), path(r2)

    output:
    path 'reads.*.fastq.gz'

    script:
    def pairs = params.query_ladder.toString().tokenize(',').join(' ')
    """
    ${params.python} ${projectDir}/mgnify_subset.py ladder --r1 ${r1} --r2 ${r2} \\
        --pairs ${pairs} --seed ${params.query_seed}
    """

    stub:
    params.query_ladder.toString().tokenize(',').collect { n -> "touch reads.${n}_1.fastq.gz reads.${n}_2.fastq.gz" }.join('\n')
}

process QUERY {
    tag "${pairs} pairs, ${name}, draws ${draws}"
    publishDir "${params.outdir}/query", mode: 'copy', pattern: '*.json'

    input:
    tuple val(pairs), path(r1), path(r2), val(name), path(index), val(draws)

    output:
    path "${name}.${pairs}.${draws}.json"

    script:
    def scratch_dir = params.query_scratch.toString() in ['', 'true', 'false'] ? '' : params.query_scratch  // CLI passes strings
    """
    export POLARS_MAX_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus}
    # Node-local copy of what the query reads: page-cache misses then cost a local read, not
    # a network one. Falls back to the index in place if the copy fails (e.g. disk full).
    idx=${index}
    if [ -n "${scratch_dir}" ] && scratch=\$(mktemp -d "${scratch_dir}/kfp-index.XXXXXX"); then
        trap 'rm -rf "\$scratch"' EXIT
        start=\$SECONDS
        if (cd -P ${index} && cp meta.json units.parquet tier2.*.npy \$(ls units.*.npy dense.*.npy 2>/dev/null) "\$scratch"); then
            idx=\$scratch
            echo "index copied to \$scratch: \$(du -sh "\$scratch" | cut -f1) in \$((SECONDS - start)) s" >&2
        else
            echo "index copy to \$scratch failed; querying ${index} in place" >&2
        fi
    fi
    if ${params.query_preload}; then  # page cache is per node, so not a task of its own
        start=\$SECONDS
        cat \$idx/tier2.*.npy > /dev/null
        echo "preload: \$(du -chL \$idx/tier2.*.npy | tail -1 | cut -f1) in \$((SECONDS - start)) s" >&2
    fi
    ${params.kfp} query \$idx ${r1} ${r2} --draws ${draws} --out profile.tsv \\
        --stats ${name}.${pairs}.${draws}.json${params.query_in_memory ? ' --in-memory' : ''}${params.query_low_memory ? ' --low-memory' : ''}
    """

    stub:
    "touch ${name}.${pairs}.${draws}.json"
}

process QUERY_COST {
    publishDir params.outdir, mode: 'copy'

    input:
    path stats

    output:
    path 'query_cost.tsv'

    script:
    """
    ${params.python} ${projectDir}/mgnify_subset.py query-cost ${stats}
    """

    stub:
    "touch query_cost.tsv"
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
    if (params.index || params.stats || params.build) {
        // Even protein_id ranges; MEMBERSHIP fails if max_protein_id leaves members out.
        def width = (params.max_protein_id as long).intdiv(params.shards as int) + 1
        MEMBERSHIP(samples[0], width)
        ch_shards = channel.of(0..<(params.shards as int)).map { i -> [i, i * width, (i + 1) * width] }
        EXTRACT(ch_shards, MEMBERSHIP.out.membership)
        MERGE(channel.of(0..<(params.buckets as int)), EXTRACT.out.parts.flatten().collect())
    }
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
        ch_part_ids = ch_parts.map { parts -> parts.collect { p -> p[0] } }
        ch_part_files = ch_parts.map { parts -> parts.collect { p -> p[1] }.flatten() }
        // Set dedup reduce: one job per set-hash range, so CONCAT holds only distinct sets.
        DEDUP(channel.of(0..<(params.ranges as int)), ch_part_ids, ch_part_files)
        CONCAT(ch_part_ids, ch_part_files, DEDUP.out.flatten().collect(), UNITS.out)
    }
    if (params.query) {
        if (!params.query_indexes) {
            error "--query needs --query_indexes, e.g. 'cost-nested/1in*/index'"
        }
        ch_run = params.query_reads
            ? channel.of(params.query_reads.toString().tokenize(',').collect { f -> file(f, checkIfExists: true) })
            : FETCH_READS(params.query_run)
        ch_reads = LADDER(ch_run).flatten()
            .map { f -> [(f.name =~ /reads\.(\d+)_/)[0][1] as long, f] }
            .groupTuple(size: 2)
            .map { n, fs -> [n] + fs.sort { f -> f.name } }
        // Index name: its parent directory (1inN).
        ch_indexes = channel.fromPath(params.query_indexes, type: 'dir', checkIfExists: true)
            .map { d -> [d.parent.name, d] }
        ch_cells = ch_reads.combine(ch_indexes)
        // draws 0 everywhere; posterior draws on one index column only.
        ch_query = ch_cells.map { c -> c + [0] }
            .mix(ch_cells.filter { c -> params.query_draws > 0 && c[3] == params.query_draws_on }
                .map { c -> c + [params.query_draws] })
        QUERY(ch_query)
        QUERY_COST(QUERY.out.collect())
    }
}
