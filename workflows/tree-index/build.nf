// Tree index, the slim, standard and large tiers (docs/tree-index-plan.md, Progress log: tiers
// by identity), built in parts. See README.md; on HPC,
// hpc/kfp-ablations/tree-index-build/run.sh. Run as: nextflow run workflows/tree-index/build.nf
//
// FIXTURE    synthetic members when --members and --slice_store are '' (the test profile)
// FAMILIES   with --slice_store: --slice_families random Pfam families of --slice_min_proteins
//            to --slice_max_proteins proteins from --slice_pfam (slice.py families)
// TAKE       per members bucket: every member of each cluster holding one of their proteins
// SLICE      the buckets' takes -> members.parquet: every unit of the families, whole
// PARTITION  proteins -> parts of ~--part_proteins (k 11 storage components and equal
//            sequences kept whole) -> parts/part*.parquet, partition.json (partition.py); with
//            --partition_per_file, once per members file (release buckets, each of whole
//            MGnify90 clusters) -> b000-part*.parquet, b000-partition.json, ...
// BUILD      per part and tier: build.py (a leaf per dense 8-mer set, dense 8-mer topology,
//            spread minimizers and the k 7 ladder diffed on the spanning tree; standard and
//            large add a protein layer on spanning-tree groups' composites, its k from the whole
//            database's size)
// MERGE      per tier: the parts -> index-<tier>/ (merge.py), query.py's format
// STORE      per part and tier: its leaves' longest proteins in leaf order (store.py) ->
//            store-<tier>/ (STORE_DIR), so verify.py and tier2.py read candidates by position
// TARGETS    --bench_targets proteins with an index position (recall.py targets)  (--bench false: skip)
// SIMULATE   one strain per target at 100-70% identity, BLOSUM62 substitutions (spike.py simulate)
// QUERY      per tier: the strains' reads, chance seeds counted on reversed reads (query.py --decoy),
//            then the candidates' rough alignment, decoy-calibrated (verify.py); tier 2: each
//            read aligned to every leaf of its routed groups, chunk by chunk (tier2.py)
// SCORE      per tier: candidate recall against the truth (recall.py score); precision and
//            read classification (recall.py reads; by Pfam with --pfam or the slice's families)
// SUMMARY    -> bench.tsv: per tier, postings, decoy bar and recall by identity (recall.py summary,
//            verify_0.01_nb: aligned at a 1% decoy bar, neighbour rule)
//
// With --reads R1,R2 and --built <dir> (a build run's outdir: index-<tier>/, store-<tier>/), a
// real sample through the built tiers instead (no build, no truth); on HPC,
// hpc/kfp-ablations/tree-index-reads/run.sh:
// TIER1      per tier: query.py --decoy, verify.py -> q-<tier>/
// TIER2      per tier and --tier2_memory budget (GB): tier2.py -> t2-<tier>-<gb>gb/ (handoff/ in
//            kfp-pileup reference's format)
// CHUNKS     per tier: the budgets' calls and assignments compared; at most 1% of reads may differ -> chunks-<tier>.json
//
// With --built <dir> alone, or --phase4 after a build, TARGETS to SCORE on built tiers and phase 4
// (docs/tree-index-plan.md, Benchmark scoring without fixed units; needs --members; on HPC,
// hpc/kfp-ablations/tree-index-benchmark/run.sh):
// FLAT       per --flat_tiers entry: today's kfp index on the same members, kfp query --aai
// REPORTED   per tree tier (tier 1 and tier 2) and flat tier: units and their sequences (nearest.py)
// ALIGN      DIAMOND: the strains against what was reported (method-independent identities)
// NEAREST    nearest-match curves, earth mover's distance, roll-ups, AAI -> nearest/<run>/
// GATE       each tree tier against the flat tier of its name -> gate.tsv

nextflow.enable.dsl = 2

params.members       = ''         // members parquet (absolute path or glob), '' = fixture
params.pfam          = ''         // protein_id, pfam_accession parquet for Pfam precision
// a family-complete slice instead of --members (docs/tree-index-plan.md, Progress log:
// precision and read classification): members buckets (glob) and the release's Pfam table
params.slice_store        = ''
params.slice_pfam         = ''
params.slice_families     = 40
params.slice_min_proteins = 20000
params.slice_max_proteins = 150000
params.part_proteins = 2000000    // proteins per part
params.partition_per_file = false // one PARTITION per members file (the release's buckets)
// a real sample through built tiers, instead of building
params.reads         = ''         // R1,R2
params.built         = ''         // a build run's outdir (index-<tier>/, store-<tier>/)
params.tier2_memory  = '16,4'     // tier2.py --memory-gb budgets; outputs must not depend on them
// decoy counts x this: bars as for a database this many times larger (query.decoy_bar)
params.chance_scale  = 1
params.store_format  = 'flat'     // the leaf-ordered store (store.py): flat (mmap) or parquet
// name:build.py options, ';'-separated (split at each tier's first ':')
params.tiers         = 'slim:--mode slim;standard:--mode standard;large:--mode large'
params.bench         = true
params.bench_targets = 2000
// sim1 and sim4 pooled (docs/tree-index-plan.md, Progress log: tiers by identity): high and low
// identities together, so the decoy bar is set on a sample with true seeds at both
params.bench_args    = '--present 1 --identities 1 0.98 0.95 0.9 0.85 0.8 0.75 0.7 --substitution blosum62 --depth-mu 0.3 --depth-sigma 0.6'

params.phase4        = false
params.flat_tiers    = 'slim:--k 11 --t-base 0.01 --t-base-singleton 0.001 --n-min 8 --max-groups 16;standard:--k 11 --t-base 0.02 --n-min 16 --oversample 1 --max-groups 16;large:--k 11 --t-base 0.05 --n-min 32 --oversample 1 --max-groups 16'
params.score_deltas  = '0.1'      // recall.py score's neighbour rule (phase 4: 0.05 0.1 0.15)
params.diamond_args  = '--more-sensitive --max-target-seqs 0 --id 50 --query-cover 50'

def py() { "${params.venv}/bin/python" }

process FIXTURE {
    label 'process_low'

    output:
    path 'members.bucket*.parquet'

    script:
    "${py()} ${projectDir}/phase0.py fixture --families 40 --members 8 --buckets 1"
}

process FAMILIES {
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy', pattern: 'families.*'

    input:
    path pfam, stageAs: 'release_pfam.parquet'
    path code, stageAs: 'code/*'

    output:
    path 'fam/chosen.parquet', emit: chosen
    path 'fam/pfam.parquet', emit: pfam
    path 'families.*'

    script:
    """
    ${py()} code/slice.py families --pfam ${pfam} --families ${params.slice_families} \
        --min-proteins ${params.slice_min_proteins} --max-proteins ${params.slice_max_proteins} \
        --threads ${task.cpus} --out fam
    cp fam/families.tsv fam/families.json .
    """
}

process TAKE {
    tag "${bucket.baseName}"
    label 'process_medium'

    input:
    tuple path(bucket), path(chosen)
    path code, stageAs: 'code/*'

    output:
    path "take.${bucket.baseName}.parquet"

    script:
    "${py()} code/slice.py take --members ${bucket} --chosen ${chosen} --threads ${task.cpus} --out take.${bucket.baseName}.parquet"
}

process SLICE {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'

    input:
    path takes, stageAs: 'takes/*'
    path code, stageAs: 'code/*'

    output:
    path 'members.parquet'

    script:
    "${py()} code/slice.py concat --parts takes/* --out members.parquet"
}

process PARTITION {
    tag "${prefix}"
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy', pattern: '*partition.json'

    input:
    tuple val(prefix), path(members, stageAs: 'members/*')
    path code, stageAs: 'code/*'
    val part_proteins  // an input, so -resume reruns when it changes

    output:
    path '*part[0-9]*.parquet', emit: parts
    path '*partition.json', emit: json

    script:
    "${py()} code/partition.py --members members/* --out . --part-proteins ${part_proteins} --prefix '${prefix}'"
}

process BUILD {
    tag "${tier} ${part.baseName}"
    label 'process_medium'

    input:
    tuple val(tier), val(args), path(part)
    val proteins  // the whole database's, for the layer's k
    path code, stageAs: 'code/*'

    output:
    tuple val(tier), path("${tier}-${part.baseName}"), emit: dir
    tuple val(tier), path("${tier}-${part.baseName}"), path(part, includeInputs: true), emit: with_part

    script:
    "${py()} code/build.py --members ${part} --out ${tier}-${part.baseName} --database-proteins ${proteins} ${args}"
}

process MERGE {
    tag "${tier}"
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', saveAs: { "index-${tier}" }

    input:
    tuple val(tier), path(parts, stageAs: 'parts/*')
    path code, stageAs: 'code/*'

    output:
    path "index-${tier}"

    script:
    "${py()} code/merge.py --parts parts/* --out index-${tier}"
}

process TARGETS {
    label 'process_medium'

    input:
    path members, stageAs: 'members/*'
    path index, stageAs: 'index'
    path code, stageAs: 'code/*'

    output:
    path 'targets.fa'

    script:
    "${py()} code/recall.py targets --members members/* --index index --n ${params.bench_targets} --out targets.fa"
}

process SIMULATE {
    label 'process_low'

    input:
    path targets
    path code, stageAs: 'code/*'

    output:
    path 'sim'

    script:
    "${py()} code/spike.py simulate --targets ${targets} ${params.bench_args} --out sim"
}

process STORE {
    tag "${tier} ${dir.name}"
    label 'process_medium'

    input:
    tuple val(tier), path(dir), path(part), path(index, stageAs: 'index')
    path code, stageAs: 'code/*'

    output:
    tuple val(tier), path('store-part/*')

    script:
    """
    ${py()} code/store.py build --proteins ${dir}/proteins.parquet --members ${part} \
        --parts index/parts.json --name ${dir.name} --format ${params.store_format} --out store-part
    """
}

process STORE_DIR {
    tag "${tier}"
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(tier), path(files, stageAs: 'files/*')

    output:
    tuple val(tier), path("store-${tier}")

    script:
    "mkdir store-${tier} && cp -L files/* store-${tier}/"
}

process QUERY {
    tag "${tier}"
    label 'process_high_memory'

    input:
    tuple val(tier), path(index, stageAs: 'index'), path(store, stageAs: 'store'), path(sim)
    path code, stageAs: 'code/*'

    output:
    tuple val(tier), path("q-${tier}")

    script:
    def reads = "${sim}/reads_1.fq.gz ${sim}/reads_2.fq.gz"
    // tier 2's chunks and batches from half the task's memory (verify.py's peak is tier 1's)
    def mem = task.memory ? task.memory.toGiga().intdiv(2) : 4
    """
    ${py()} code/query.py --index index --reads ${reads} --decoy --chance-scale ${params.chance_scale} --out q-${tier}
    ${py()} code/verify.py --index index --store store --reads ${reads} --query q-${tier} --chance-scale ${params.chance_scale}
    ${py()} code/tier2.py --index index --store store --reads ${reads} --query q-${tier} --chance-scale ${params.chance_scale} \
        --memory-gb ${mem} --out q-${tier}/tier2
    """
}

process SCORE {
    tag "${tier}"
    label 'process_medium'

    input:
    tuple val(tier), path(index, stageAs: 'index'), path(q), path(sim)
    path pfam, stageAs: 'pfam/*'
    path code, stageAs: 'code/*'

    output:
    tuple val(tier), path("s-${tier}")

    script:
    def fam = pfam ? "--pfam ${pfam[0]}" : ''
    """
    ${py()} code/recall.py score --index index --query ${q} --truth ${sim}/truth.tsv \
        --out s-${tier} --caps 256 --deltas ${params.score_deltas}
    ${py()} code/recall.py reads --index index --query ${q} --truth ${sim}/truth.tsv \
        --read-truth ${sim}/reads.parquet ${fam} --out s-${tier}
    ${py()} code/recall.py reads --index index --query ${q}/tier2 --truth ${sim}/truth.tsv \
        --read-truth ${sim}/reads.parquet ${fam} --out s-${tier}/tier2
    """
}

process TIER1 {
    tag "${tier}"
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy', pattern: 'q-*'

    input:
    tuple val(tier), path(index, stageAs: 'index'), path(store, stageAs: 'store')
    path reads, stageAs: 'reads/*'
    path code, stageAs: 'code/*'

    output:
    tuple val(tier), path(index), path(store), path("q-${tier}")

    script:
    """
    ${py()} code/query.py --index index --reads ${reads} --decoy --chance-scale ${params.chance_scale} --out q-${tier}
    ${py()} code/verify.py --index index --store store --reads ${reads} --query q-${tier} --chance-scale ${params.chance_scale}
    """
}

process TIER2 {
    tag "${tier} ${gb} GB"
    label 'process_high_memory'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(tier), path(index, stageAs: 'index'), path(store, stageAs: 'store'), path(q), val(gb)
    path reads, stageAs: 'reads/*'
    path code, stageAs: 'code/*'

    output:
    tuple val(tier), path("t2-${tier}-${gb}gb")

    script:
    """
    ${py()} code/tier2.py --index index --store store --reads ${reads} --query ${q} --chance-scale ${params.chance_scale} \
        --memory-gb ${gb} --out t2-${tier}-${gb}gb
    rm -rf t2-${tier}-${gb}gb/spill
    """
}

process CHUNKS {
    tag "${tier}"
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(tier), path(runs, stageAs: 'runs/*')
    path code, stageAs: 'code/*'

    output:
    path "chunks-${tier}.json"

    script:
    "${py()} code/recall.py chunks --runs runs/* --out chunks-${tier}.json"
}

process SUMMARY {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path runs, stageAs: 'runs/*'
    path code, stageAs: 'code/*'

    output:
    path 'bench.tsv'
    path 'runs/[qs]-*', includeInputs: true

    script:
    "${py()} code/recall.py summary --root runs --base base --strategy verify_0.01_nb --out bench.tsv"
}

process FLAT {
    tag "${tier}"
    label 'process_high_memory'

    input:
    tuple val(tier), val(index_args)
    path members, stageAs: 'members/*'
    path sim

    output:
    tuple val("flat-${tier}"), path("flat-${tier}")

    script:
    """
    mkdir flat-${tier}
    ${params.kfp} index members/* flat-${tier}/index ${index_args}
    ${params.kfp} query flat-${tier}/index ${sim}/reads_1.fq.gz ${sim}/reads_2.fq.gz --aai \
        --out flat-${tier}/profile.tsv --stats flat-${tier}/stats.json
    """
}

process REPORTED {
    tag "${run}"
    label 'process_medium'

    input:
    tuple val(run), val(args), path(data, stageAs: 'data/*'), path(sim)
    path members, stageAs: 'members/*'
    path code, stageAs: 'code/*'

    output:
    tuple val(run), path("rep-${run}")

    script:
    "${py()} code/nearest.py reported --truth ${sim}/truth.tsv --members members/* ${args} --out rep-${run}"
}

process ALIGN {
    tag "${run}"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.2.8--he361c42_0'

    input:
    tuple val(run), path(rep)

    output:
    tuple val(run), path('hits.tsv')

    script:
    """
    diamond makedb --in ${rep}/reported.fa --db reported --threads ${task.cpus}
    diamond blastp --db reported --query ${rep}/strains.fa --out hits.tsv --threads ${task.cpus} \
        --outfmt 6 qseqid sseqid pident ${params.diamond_args}
    """
}

process NEAREST {
    tag "${run}"
    label 'process_medium'
    publishDir "${params.outdir}/nearest", mode: 'copy'

    input:
    tuple val(run), path(rep), path(hits), path(extra, stageAs: 'extra/*'), path(sim)
    path members, stageAs: 'members/*'
    path pfam, stageAs: 'pfam/*'
    path code, stageAs: 'code/*'

    output:
    path run

    script:
    def fam = pfam ? "--rollup pfam=pfam/*" : ''
    def strains = extra ? "--strains extra/*/strains.parquet" : ''
    """
    ${py()} code/nearest.py score --reported ${rep} --truth ${sim}/truth.tsv --hits ${hits} \
        --members members/* ${fam} ${strains} --out ${run}
    """
}

process GATE {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    input:
    path runs
    path code, stageAs: 'code/*'

    output:
    path 'gate.tsv'

    script:
    "${py()} code/nearest.py gate ${runs} --out gate.tsv"
}

workflow {
    code = ['partition.py', 'build.py', 'merge.py', 'survey.py', 'phase0.py', 'query.py',
            'alphabets.py', 'recall.py', 'verify.py', 'tier2.py', 'store.py', 'report.py', 'slice.py',
            '../spike-benchmark/spike.py']
        .collect { file("${projectDir}/${it}") }
    pfam = params.pfam ? Channel.fromPath(params.pfam, checkIfExists: true).collect()
                       : Channel.value([])
    tiers = Channel.fromList(params.tiers.tokenize(';')*.trim().collect { it.split(':', 2) as List })
    if (params.reads) {  // a real sample through built tiers
        reads = params.reads.tokenize(',').collect { file(it, checkIfExists: true) }
        built = tiers.map { tier, args ->
            [tier, file("${params.built}/index-${tier}", checkIfExists: true),
             file("${params.built}/store-${tier}", checkIfExists: true)] }
        TIER1(built, reads, code)
        TIER2(TIER1.out.combine(Channel.fromList(params.tier2_memory.toString().tokenize(','))), reads, code)
        CHUNKS(TIER2.out.groupTuple(), code)
    } else if (params.built) {  // phase 4 on built tiers
        members = Channel.fromPath(params.members, checkIfExists: true).collect()
        index = tiers.map { tier, args -> [tier, file("${params.built}/index-${tier}", checkIfExists: true)] }
        stores = tiers.map { tier, args -> [tier, file("${params.built}/store-${tier}", checkIfExists: true)] }
        bench(index, stores, members, pfam, code, true)
    } else {
        build(code, pfam, tiers)
    }
}

workflow build {
    take:
    code
    pfam
    tiers

    main:
    if (params.slice_store) {
        FAMILIES(Channel.fromPath(params.slice_pfam, checkIfExists: true).first(), code)
        TAKE(Channel.fromPath(params.slice_store, checkIfExists: true).combine(FAMILIES.out.chosen), code)
        members = SLICE(TAKE.out.collect(), code).collect()
        pfam = FAMILIES.out.pfam.collect()
    } else {
        members = params.members ? Channel.fromPath(params.members, checkIfExists: true).collect()
                                 : FIXTURE()
    }
    // one partition, or one per members file (release buckets: b000-, b001-, ...)
    by_file = members.flatMap { fs ->
        fs = fs instanceof List ? fs : [fs]
        params.partition_per_file
            ? fs.sort { it.name }.withIndex().collect { f, i -> [String.format('b%03d-', i), [f]] }
            : [['', fs]]
    }
    PARTITION(by_file, code, params.part_proteins)
    proteins = PARTITION.out.json.map { new groovy.json.JsonSlurper().parse(it).proteins as long }.sum()
    BUILD(tiers.combine(PARTITION.out.parts.flatten()), proteins, code)
    MERGE(BUILD.out.dir.groupTuple(), code)
    index = MERGE.out.map { [it.name - 'index-', it] }
    // the leaf-ordered store, per part (its positions one range of the merged index)
    STORE(BUILD.out.with_part.combine(index, by: 0), code)
    stores = STORE_DIR(STORE.out.groupTuple().map { [it[0], it[1].flatten()] })
    if (params.bench) {
        bench(index, stores, members, pfam, code, params.phase4)
    }
}

workflow bench {
    take:
    index
    stores
    members
    pfam
    code
    phase4

    main:
    TARGETS(members, index.first().map { it[1] }, code)
    sim = SIMULATE(TARGETS.out, code).first()
    QUERY(index.join(stores).combine(sim), code)
    SCORE(index.join(QUERY.out).combine(sim), pfam, code)
    // index-<tier>, q-<tier>, s-<tier> side by side
    SUMMARY(index.map { it[1] }.mix(QUERY.out.map { it[1] }, SCORE.out.map { it[1] }).collect(), code)
    if (phase4) {
        scorer = code + [file("${projectDir}/nearest.py")]
        FLAT(Channel.fromList(params.flat_tiers.tokenize(';')*.trim().collect { it.split(':', 2) as List }),
             members, sim)
        tree = index.join(QUERY.out).flatMap { tier, i, q ->
            def idx = "--index data/${i.name} --index-dir data/${i.name}"
            [["tree-${tier}", "${idx} --query data/${q.name}", [i, q]],
             ["tree-${tier}-tier2", "${idx} --query data/${q.name}/tier2", [i, q]]] }
        flat = FLAT.out.map { run, d ->
            [run, "--profile data/${d.name}/profile.tsv --stats data/${d.name}/stats.json --index-dir data/${d.name}/index", [d]] }
        REPORTED(tree.mix(flat).combine(sim), members, scorer)
        ALIGN(REPORTED.out)
        // recall.py score's strains.parquet for the tree's runs, for the neighbour rule
        strains = REPORTED.out.join(ALIGN.out).combine(SCORE.out).filter { run, rep, hits, tier, s ->
            run == "tree-${tier}" || run == "tree-${tier}-tier2" }.map { run, rep, hits, tier, s -> [run, rep, hits, [s]] }
        others = REPORTED.out.join(ALIGN.out).filter { it[0].startsWith('flat-') }.map { it + [[]] }
        NEAREST(strains.mix(others).combine(sim), members, pfam, scorer)
        GATE(NEAREST.out.collect(), scorer)
    }
}
