// fmh-funprofiler benchmark: KO and Pfam detection by kmer-functional-profiler indexes and
// by the fmh-funprofiler KO sketches (compat mode) on InSilicoSeq metagenomes, and KO
// detection by other tools (--tools: DIAMOND, fmh-funprofiler, kMermaid, HUMAnN 3 and 4).
// See README.md.


process FETCH {
    label 'process_single'
    storeDir params.data_dir

    output:
    path 'genomes_extracted_from_kegg', emit: genomes
    path 'protein_ref_db_giant.faa',    emit: faa
    path 'present_genes_and_koids.csv', emit: kos
    path "${params.sketches}",          emit: sketches

    script:
    """
    Z=https://zenodo.org/records/10055954/files
    curl -fsSL -o genomes.zip "\$Z/genomes_extracted_from_kegg.zip?download=1"
    unzip -q genomes.zip -x '*.DS_Store' '__MACOSX/*' && rm genomes.zip
    for f in protein_ref_db_giant.faa present_genes_and_koids.csv ${params.sketches}; do
        curl -fsSL -o "\$f" "\$Z/\$f?download=1"
    done
    """

    stub:
    "mkdir genomes_extracted_from_kegg && touch protein_ref_db_giant.faa present_genes_and_koids.csv ${params.sketches}"
}

process MEMBERS {
    label 'process_medium'

    input:
    path faa
    path kos

    output:
    path 'members.parquet', emit: members

    script:
    "${params.bench} members --faa ${faa} --kos ${kos}"

    stub:
    "touch members.parquet"
}

process INDEX {
    tag "${name}"
    label 'process_high_memory'
    // meta.json is copied out of the index: a file inside a declared directory output is
    // folded into that directory and never matches the publish pattern
    publishDir params.outdir, mode: 'copy', pattern: 'meta.json', saveAs: { "index_${name}/meta.json" }

    input:
    tuple val(label), val(name), val(args), path(members)
    path code, stageAs: 'code/*'  // package sources: only here so -resume rebuilds on changes

    output:
    tuple val(label), val(name), path("index_${name}"), emit: index
    path 'meta.json'

    script:
    "${params.kfp} index ${members} index_${name} ${args} && cp index_${name}/meta.json meta.json"

    stub:
    "mkdir index_${name} && touch index_${name}/meta.json meta.json"
}

process IMPORT_SKETCHES {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', pattern: 'meta.json', saveAs: { 'index_fmh_compat/meta.json' }

    input:
    path sketches
    path code, stageAs: 'code/*'  // package sources: only here so -resume rebuilds on changes

    output:
    tuple val('ko'), val('fmh_compat'), path('index_fmh_compat'), emit: index
    path 'meta.json'

    script:
    """
    ${params.kfp} import-sourmash ${sketches} index_fmh_compat --ksize ${params.ksize}
    cp index_fmh_compat/meta.json meta.json
    """

    stub:
    "mkdir index_fmh_compat && touch index_fmh_compat/meta.json meta.json"
}

// ---- Pfam labels (Benchmark labels in the plan): the genomes' proteins annotated with
// Pfam-A (hmmsearch --cut_ga), domains as truth features and as members of Pfam units.

process PFAM_DB {
    label 'process_single'
    storeDir "${params.db_dir}/pfam"

    output:
    path 'Pfam-A.hmm', emit: hmm
    path 'Pfam.version'

    script:
    """
    curl -fsSL ${params.pfam_url}/Pfam-A.hmm.gz | gunzip > Pfam-A.hmm
    curl -fsSL ${params.pfam_url}/Pfam.version.gz | gunzip > Pfam.version
    """

    stub:
    "touch Pfam-A.hmm Pfam.version"
}

process PFAM_PROTEINS {
    label 'process_single'

    input:
    path genomes
    path code, stageAs: 'code/*'  // bench.py: only here so -resume reruns on changes

    output:
    path 'proteins.*.faa', emit: faa

    script:
    "${params.bench} pfam-proteins --genomes-dir ${genomes} --chunks ${params.pfam_chunks}"

    stub:
    "touch proteins.0.faa proteins.1.faa"
}

process PFAM_ANNOTATE {
    tag "${faa.baseName}"
    label 'process_medium'
    container 'quay.io/biocontainers/hmmer:3.4--hdbdd923_2'

    input:
    path faa
    path hmm

    output:
    path "${faa.baseName}.domtbl", emit: domtbl

    script:
    """
    hmmsearch --cpu ${task.cpus} ${params.pfam_threshold} --domtblout ${faa.baseName}.domtbl \\
        -o /dev/null ${hmm} ${faa}
    """

    stub:
    "touch ${faa.baseName}.domtbl"
}

process PFAM_DOMAINS {
    label 'process_medium'
    publishDir "${params.outdir}/pfam", mode: 'copy'

    input:
    path domtbl, stageAs: 'tbl/*'
    path faa, stageAs: 'faa/*'
    path code, stageAs: 'code/*'  // bench.py: only here so -resume reruns on changes

    output:
    path 'domains.parquet', emit: domains
    path 'pfam_members.parquet', emit: members

    script:
    "${params.bench} pfam-domains --domtbl tbl/* --proteins faa/*"

    stub:
    "touch domains.parquet pfam_members.parquet"
}

process SAMPLE {
    tag "seed ${seed}"
    label 'process_single'

    input:
    val seed
    path genomes

    output:
    tuple val(seed), path('sample.fna'), path('genes.parquet'), emit: sample
    path 'genomes.txt'

    script:
    "${params.bench} sample --genomes-dir ${genomes} --n ${params.n_genomes} --seed ${seed}"

    stub:
    "touch sample.fna genes.parquet genomes.txt"
}

process SIMULATE {
    tag "seed ${seed}"
    label 'process_medium'

    input:
    tuple val(seed), path(fna)

    output:
    tuple val(seed), path('reads_R1.fastq.gz'), path('reads_R2.fastq.gz'), emit: reads

    script:
    // bench.py iss patches the perfect model's bugs in iss 2.0.1
    """
    ${params.iss_mode == 'perfect' ? "${params.bench} iss generate --mode perfect" : "${params.venv}/bin/iss generate --model ${params.iss_model}"} --genomes ${fna} \\
        --n_reads ${params.n_reads} --abundance lognormal --seed ${seed} \\
        --cpus ${task.cpus} --compress --output reads
    """

    stub:
    "touch reads_R1.fastq.gz reads_R2.fastq.gz"
}

process TRUTH {
    tag "seed ${seed}"
    label 'process_medium'
    publishDir "${params.outdir}/truth", mode: 'copy', saveAs: { f -> f == 'truth.csv' ? "seed${seed}.csv" : "seed${seed}_pfam.csv" }

    input:
    tuple val(seed), path(fna), path(genes), path(r1), path(r2)
    path kos
    path domains  // [] without Pfam
    path code, stageAs: 'code/*'  // bench.py: only here so -resume reruns on changes

    output:
    tuple val(seed), val('ko'), path('truth.csv'), emit: truth
    tuple val(seed), val('pfam'), path('truth_pfam.csv'), emit: pfam, optional: true

    script:
    """
    ${params.bench} truth --fna ${fna} --genes ${genes} --kos ${kos} --r1 ${r1} --r2 ${r2} \\
        --threads ${task.cpus} ${domains ? "--domains ${domains}" : ''}
    """

    stub:
    "touch truth.csv ${domains ? 'truth_pfam.csv' : ''}"
}

process PROFILE {
    tag "seed ${seed} ${name}"
    label 'process_single'
    // one saveAs for both outputs: without the per-file branch they overwrite each other
    publishDir params.outdir, mode: 'copy', saveAs: { f ->
        f.endsWith('.parquet') ? "kmers/seed${seed}_${name}.parquet" : "profiles/seed${seed}_${name}.tsv"
    }

    input:
    tuple val(seed), path(r1), path(r2), val(label), val(name), path(index), val(by_pfam)
    path code, stageAs: 'code/*'  // query sources: only here so -resume reruns on changes

    output:
    tuple val(seed), val(label), val(name), path('profile.tsv'), path('kmers.parquet'), emit: profile

    script:
    // by_pfam: units carry Pfam labels (MGnify90 clusters), and the profile is summed per Pfam
    def out = by_pfam ? 'units.tsv' : 'profile.tsv'
    def per_pfam = "${params.bench} pfam-profile --profile units.tsv --unit-pfam ${index}/unit_pfam.parquet"
    """
    ${params.kfp} query ${index} ${r1} ${r2} --out ${out} --draws ${params.draws} --kmers kmers.parquet --all-estimators
    ${by_pfam ? "${per_pfam} --out profile.tsv" : ''}
    """

    stub:
    "touch profile.tsv kmers.parquet"
}

process SCORE {
    tag "seed ${seed} ${name}"
    label 'process_single'

    input:
    tuple val(seed), val(label), val(name), path(profile), path(truth)
    path code, stageAs: 'code/*'  // bench.py: only here so -resume reruns on changes

    output:
    path 'score.tsv', emit: score

    script:
    // DIAMOND's evidence is read pairs, not k-mers: its own range, to trace a curve
    def min_hits = name == 'diamond' ? params.diamond_min_hits : params.min_hits
    """
    ${params.bench} score --truth ${truth} --profile ${profile} --sample seed${seed} --index ${name} \\
        --label ${label} --min-hits ${min_hits.toString().tokenize(',').join(' ')}
    """

    stub:
    "touch score.tsv"
}

process DETECTED {
    tag "seed ${seed} ${name}"
    label 'process_single'
    publishDir "${params.outdir}/detected", mode: 'copy', saveAs: { "seed${seed}_${name}.tsv" }

    input:
    tuple val(label), val(name), val(seed), path(profile), path(kmers), path(truth), path(fna), path(index)
    path code, stageAs: 'code/*'  // bench.py: only here so -resume reruns on changes

    output:
    path 'detected.tsv', emit: detected

    script:
    """
    ${params.bench} detected --truth ${truth} --profile ${profile} --kmers ${kmers} \\
        --index-dir ${index} --fna ${fna} --sample seed${seed} --index ${name}
    """

    stub:
    "touch detected.tsv"
}

// ---- Other tools (--tools). Each writes its raw output to raw/; TOOL_PROFILE turns it into
// a profile SCORE reads. Tools run in pinned containers
// (-profile docker or singularity); the kMermaid and HUMAnN images are built by
// containers/build.sh.

process DIAMOND_DB {
    label 'process_medium'
    storeDir "${params.db_dir}/diamond"
    container 'quay.io/biocontainers/diamond:2.2.8--he361c42_0'

    input:
    path faa

    output:
    path 'kegg.dmnd', emit: db

    script:
    "diamond makedb --in ${faa} --db kegg --threads ${task.cpus}"

    stub:
    "touch kegg.dmnd"
}

process DIAMOND {
    tag "seed ${seed}"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.2.8--he361c42_0'

    input:
    tuple val(seed), path(r1), path(r2)
    path db

    output:
    tuple val(seed), val('diamond'), path('raw'), emit: raw

    script:
    // best hit per read (-k 1) against the KEGG proteins the KO indexes are built from
    """
    mkdir raw
    cat ${r1} ${r2} > reads.fastq.gz
    diamond blastx --db ${db} --query reads.fastq.gz --out raw/hits.tsv --threads ${task.cpus} \\
        --max-target-seqs 1 --outfmt 6 qseqid sseqid length slen bitscore ${params.diamond_args}
    rm reads.fastq.gz
    """

    stub:
    "mkdir raw && touch raw/hits.tsv"
}

process FMH_FUNPROFILER {
    tag "seed ${seed}"
    label 'process_single'
    container 'quay.io/biocontainers/fmh-funprofiler:1.1.1--pyh106432d_0'

    input:
    tuple val(seed), path(r1), path(r2)
    path sketches

    output:
    tuple val(seed), val('fmh_funprofiler'), path('raw'), emit: raw

    script:
    // the released tool as its README runs it, on one FASTQ per sample. --threshold_bp =
    // scaled is one shared hash, its default (1000) at scaled 1000
    """
    mkdir raw
    cat ${r1} ${r2} > reads.fastq.gz
    funcprofiler reads.fastq.gz ${sketches} ${params.ksize} ${params.sketch_scaled} raw/ko.csv \\
        -p raw/prefetch.csv -t ${params.sketch_scaled}
    rm reads.fastq.gz reads.fastq.gz_sketch_*.sig.zip
    """

    stub:
    "mkdir raw && touch raw/ko.csv raw/prefetch.csv"
}

process KMERMAID_MODEL {
    label 'process_high_memory'
    storeDir "${params.db_dir}/kmermaid_max${params.kmermaid_max_members}"
    container params.kmermaid_container

    input:
    path members

    output:
    path 'kmermaid_model', emit: model

    script:
    """
    kmermaid_kfp.py train --members ${members} \\
        --max-members ${params.kmermaid_max_members} > train.log
    """

    stub:
    "mkdir kmermaid_model"
}

process KMERMAID {
    tag "seed ${seed}"
    label 'process_kmermaid'
    container params.kmermaid_container

    input:
    tuple val(seed), path(r1), path(r2)
    path model

    output:
    tuple val(seed), val('kmermaid'), path('raw'), emit: raw

    script:
    // single-threaded Python; mates are classified separately, as kMermaid takes one file
    """
    mkdir raw
    zcat ${r1} ${r2} > reads.fastq
    kmermaid_kfp.py classify --model ${model} --reads reads.fastq \\
        --out raw/kmermaid.tsv
    rm reads.fastq
    """

    stub:
    "mkdir raw && touch raw/kmermaid.tsv"
}

process HUMANN_DB {
    tag "${tool} ${db}"
    label 'process_medium'
    storeDir "${params.db_dir}/${tool}"  // one database per task: a failed download loses only it
    container "${image}"

    input:
    tuple val(tool), val(image), val(db), val(build)

    output:
    tuple val(tool), val(db), path(db), emit: db

    script:
    "humann_databases --download ${db} ${build} . --update-config no"

    stub:
    "mkdir ${db}"
}

process METAPHLAN_DB {
    tag "${tool}"
    label 'process_medium'
    storeDir "${params.db_dir}/${tool}"
    container "${image}"

    input:
    tuple val(tool), val(image), val(index)

    output:
    tuple val(tool), path('metaphlan'), emit: db

    script:
    """
    metaphlan --install --index ${index} --bowtie2db metaphlan --nproc ${task.cpus}
    """

    stub:
    "mkdir metaphlan"
}

process HUMANN {
    tag "seed ${seed}"
    label 'process_humann'
    container params.humann3_container  // MetaPhlAn 4.1.1, DIAMOND 2.0.15, bowtie2 2.5.5

    input:
    tuple val(seed), path(r1), path(r2)
    tuple path(chocophlan), path(uniref), path(utility), path(mpa)

    output:
    tuple val(seed), val('humann'), path('raw'), emit: raw

    script:
    // HUMAnN takes unpaired reads: both mates in one file, as its docs advise. UniRef90
    // families are regrouped to KOs with HUMAnN's own UniRef90 -> KO mapping.
    """
    mkdir raw
    cat ${r1} ${r2} > reads.fastq.gz
    humann --input reads.fastq.gz --output out --threads ${task.cpus} \\
        --nucleotide-database ${chocophlan} --protein-database ${uniref} \\
        --metaphlan-options "--bowtie2db ${mpa} --index ${params.metaphlan_index} --nproc ${task.cpus}"
    humann_regroup_table --input out/reads_genefamilies.tsv \\
        --custom ${utility}/map_ko_uniref90.txt.gz --output raw/ko.tsv
    cp out/reads_genefamilies.tsv raw/
    rm -r reads.fastq.gz out/reads_humann_temp
    """

    stub:
    "mkdir raw && touch raw/ko.tsv"
}

process HUMANN4 {
    tag "seed ${seed}"
    label 'process_humann'
    container params.humann4_container

    input:
    tuple val(seed), path(r1), path(r2)
    tuple path(chocophlan), path(uniref), path(utility), path(mpa)

    output:
    tuple val(seed), val('humann4'), path('raw'), emit: raw

    script:
    // MetaPhlAn runs first, as HUMAnN would (same -t), and its profile goes in with
    // --taxonomic-profile: HUMAnN 4 alpha.2 otherwise looks for its database tag in
    // `metaphlan --version`, which MetaPhlAn 4.1.1 does not print, and exits.
    """
    mkdir raw
    cat ${r1} ${r2} > reads.fastq.gz
    metaphlan reads.fastq.gz --input_type fastq --bowtie2db ${mpa} \\
        --index ${params.humann4_metaphlan_index} --nproc ${task.cpus} -t rel_ab_w_read_stats \\
        --bowtie2out metaphlan.bowtie2.bz2 -o raw/metaphlan.tsv
    humann --input reads.fastq.gz --output out --threads ${task.cpus} \\
        --taxonomic-profile raw/metaphlan.tsv \\
        --nucleotide-database ${chocophlan} --protein-database ${uniref} \\
        --utility-database ${utility}
    humann_regroup_table --input out/reads_2_genefamilies.tsv \\
        --custom ${utility}/map_ko_uniref90.txt.gz --output raw/ko.tsv
    cp out/reads_2_genefamilies.tsv raw/
    rm -r reads.fastq.gz metaphlan.bowtie2.bz2 out/reads_humann_temp
    """

    stub:
    "mkdir raw && touch raw/ko.tsv"
}

process TOOL_PROFILE {
    tag "seed ${seed} ${tool}"
    label 'process_single'
    publishDir "${params.outdir}/profiles", mode: 'copy', saveAs: { "seed${seed}_${tool}.tsv" }

    input:
    tuple val(seed), val(tool), path(raw)
    path kos
    path members
    path code, stageAs: 'code/*'  // bench.py: only here so -resume reruns on changes

    output:
    tuple val(seed), val(tool), path('profile.tsv'), emit: profile

    script:
    """
    ${params.bench} tool-profile --tool ${tool} --raw ${raw} --kos ${kos} --members ${members} \\
        --scaled ${params.sketch_scaled}
    """

    stub:
    "touch profile.tsv"
}

process SUMMARY {
    label 'process_single'
    publishDir params.outdir, mode: 'copy'

    input:
    path scores, stageAs: 'score*.tsv'

    output:
    path 'summary.tsv'
    path 'scores.tsv'

    script:
    "${params.bench} summary ${scores}"

    stub:
    "touch summary.tsv scores.tsv"
}

// Output of a git command in the pipeline's checkout (for run.json).
def git(List cmd) {
    return (['git', '-C', projectDir.toString()] + cmd).execute().text.trim()
}

workflow {
    FETCH()
    MEMBERS(FETCH.out.faa, FETCH.out.kos)
    def bench_py = file("${projectDir}/bench.py")
    // ponytail: tracks the Python package only; Rust kernel changes still need a fresh run
    ch_code = channel.fromPath("${projectDir}/../../python/kmer_functional_profiler/*.py").collect()
    def labels = params.labels.toString().tokenize(',')*.trim()
    if (labels - ['ko', 'pfam']) {
        error "Unknown --labels: ${(labels - ['ko', 'pfam']).join(', ')}"
    }
    def ch_members = channel.of('ko').combine(MEMBERS.out.members)
    def ch_domains = channel.value([])
    if ('pfam' in labels) {
        def ch_hmm = params.pfam_hmm ? channel.value(file(params.pfam_hmm)) : PFAM_DB().hmm
        PFAM_PROTEINS(FETCH.out.genomes, bench_py)
        PFAM_ANNOTATE(PFAM_PROTEINS.out.faa.flatten(), ch_hmm)
        PFAM_DOMAINS(PFAM_ANNOTATE.out.domtbl.collect(), PFAM_PROTEINS.out.faa, bench_py)
        ch_domains = PFAM_DOMAINS.out.domains
        ch_members = ch_members.mix(channel.of('pfam').combine(PFAM_DOMAINS.out.members))
    }
    // Every index config is built per label: KO units (names as before) and Pfam units
    // (pfam_<name>), each Pfam domain of the genomes' proteins a member of its Pfam's unit
    ch_configs = channel.fromList(params.indexes)
        .combine(ch_members)
        .filter { _cfg, label, _members -> label in labels }
        .map { cfg, label, members -> [label, label == 'ko' ? cfg.name : "${label}_${cfg.name}", cfg.args, members] }
    INDEX(ch_configs, ch_code)
    IMPORT_SKETCHES(FETCH.out.sketches, ch_code)
    // label, name, index, by_pfam (the profile is summed per Pfam through unit_pfam.parquet)
    ch_indexes = INDEX.out.index.mix(IMPORT_SKETCHES.out.index)
        .filter { it[0] in labels }
        .map { it + [false] }
        .mix(channel.fromList(params.mgnify_indexes).filter { 'pfam' in labels }
            .map { cfg -> ['pfam', cfg.name, file(cfg.path, checkIfExists: true), true] })

    SAMPLE(channel.of(1..params.replicates), FETCH.out.genomes)
    SIMULATE(SAMPLE.out.sample.map { seed, fna, _genes -> [seed, fna] })
    TRUTH(SAMPLE.out.sample.join(SIMULATE.out.reads), FETCH.out.kos, ch_domains, bench_py)
    ch_truth = TRUTH.out.truth.mix(TRUTH.out.pfam).filter { it[1] in labels }  // seed, label, truth

    PROFILE(SIMULATE.out.reads.combine(ch_indexes), ch_code)
    ch_scored = PROFILE.out.profile.combine(ch_truth, by: [0, 1])  // seed, label, name, profile, kmers, truth
    DETECTED(
        ch_scored
            .combine(SAMPLE.out.sample.map { seed, fna, _genes -> [seed, fna] }, by: 0)
            .map { seed, label, name, profile, kmers, truth, fna -> [label, name, seed, profile, kmers, truth, fna] }
            .combine(ch_indexes.filter { !it[3] }.map { it.take(3) }, by: [0, 1]),
        bench_py,
    )
    DETECTED.out.detected.collectFile(name: 'detected.tsv', keepHeader: true, storeDir: params.outdir)

    // `--tools ''` arrives as true: no tools
    def tools = params.tools.toString() in ['', 'true', 'false'] ? [] : params.tools.toString().tokenize(',')*.trim()
    def unknown = tools - ['diamond', 'fmh_funprofiler', 'kmermaid', 'humann', 'humann4']
    if (unknown) {
        error "Unknown --tools: ${unknown.join(', ')}"
    }
    def ch_raw = channel.empty()
    if ('diamond' in tools) {
        DIAMOND(SIMULATE.out.reads, DIAMOND_DB(FETCH.out.faa).db)
        ch_raw = ch_raw.mix(DIAMOND.out.raw)
    }
    if ('fmh_funprofiler' in tools) {
        FMH_FUNPROFILER(SIMULATE.out.reads, FETCH.out.sketches)
        ch_raw = ch_raw.mix(FMH_FUNPROFILER.out.raw)
    }
    if ('kmermaid' in tools) {
        KMERMAID(SIMULATE.out.reads, KMERMAID_MODEL(MEMBERS.out.members).model)
        ch_raw = ch_raw.mix(KMERMAID.out.raw)
    }
    // HUMAnN 3.9 (~45 GB) and 4 alpha (~70 GB): ChocoPhlAn full, a UniRef90 DIAMOND database
    // (4 distributes only the EC-filtered one), utility mapping, and the MetaPhlAn index each uses
    def humann = [
        humann:  [params.humann3_container, 'uniref90_diamond', params.metaphlan_index],
        humann4: [params.humann4_container, 'uniref90_ec_filtered_diamond', params.humann4_metaphlan_index],
    ].findAll { tool, _cfg -> tool in tools }
    def dbs = ['chocophlan', 'uniref', 'utility_mapping']
    HUMANN_DB(channel.fromList(humann.collectMany { tool, cfg ->
        [[tool, cfg[0], 'chocophlan', 'full'], [tool, cfg[0], 'uniref', cfg[1]], [tool, cfg[0], 'utility_mapping', 'full']]
    }))
    METAPHLAN_DB(channel.fromList(humann.collect { tool, cfg -> [tool, cfg[0], cfg[2]] }))
    // tool -> [chocophlan, uniref, utility_mapping, metaphlan], as HUMANN and HUMANN4 take them
    def ch_humann_db = HUMANN_DB.out.db
        .groupTuple(size: dbs.size())
        .join(METAPHLAN_DB.out.db)
        .map { tool, names, paths, mpa -> [tool] + dbs.collect { db -> paths[names.indexOf(db)] } + [mpa] }
    if ('humann' in tools) {
        HUMANN(SIMULATE.out.reads, ch_humann_db.filter { it[0] == 'humann' }.map { it.drop(1) }.first())
        ch_raw = ch_raw.mix(HUMANN.out.raw)
    }
    if ('humann4' in tools) {
        HUMANN4(SIMULATE.out.reads, ch_humann_db.filter { it[0] == 'humann4' }.map { it.drop(1) }.first())
        ch_raw = ch_raw.mix(HUMANN4.out.raw)
    }
    TOOL_PROFILE(ch_raw, FETCH.out.kos, MEMBERS.out.members, bench_py)

    SCORE(
        ch_scored
            .map { seed, label, name, profile, _kmers, truth -> [seed, label, name, profile, truth] }
            .mix(TOOL_PROFILE.out.profile.map { seed, tool, profile -> [seed, 'ko', tool, profile] }
                .combine(TRUTH.out.truth, by: [0, 1])),
        bench_py,
    )
    SUMMARY(SCORE.out.score.collect())

    // run.json: what produced the results in outdir (reads, indexes, code, status). params and
    // workflow are read here: inside the handler, names resolve against the workflow metadata.
    def run_params = params
    def wf = workflow
    def reads = params.iss_mode == 'perfect' ? 'error-free (iss perfect)' : "iss ${params.iss_model} (${params.iss_mode})"
    def description = "${params.replicates} metagenomes x ${params.n_genomes} KEGG genomes, " +
        "${params.n_reads} reads, ${reads}; labels ${params.labels}; ${params.indexes.size()} index configs" +
        " + fmh_compat + ${params.mgnify_indexes.size()} MGnify; draws ${params.draws}"
    def code = [commit: git(['rev-parse', 'HEAD']), branch: git(['rev-parse', '--abbrev-ref', 'HEAD']),
                uncommitted_changes: !git(['status', '--porcelain']).isEmpty()]
    def out = file(params.outdir)
    wf.onComplete {
        def info = [
            description: description,
            success: wf.success,
            exit_status: wf.exitStatus,
            start: wf.start.toString(),
            complete: wf.complete.toString(),
            duration: wf.duration.toString(),
            command_line: wf.commandLine,
            profile: wf.profile,
            resume: wf.resume,
            session_id: wf.sessionId.toString(),
            run_name: wf.runName,
            nextflow: wf.nextflow.version.toString(),
            code: code,
            params: run_params,
        ]
        out.mkdirs()
        out.resolve('run.json').text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info)) + '\n'
    }
}
