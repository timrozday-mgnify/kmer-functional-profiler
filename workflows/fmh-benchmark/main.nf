// fmh-funprofiler benchmark: KO and Pfam detection by kmer-functional-profiler indexes and
// by the fmh-funprofiler KO sketches (compat mode) on InSilicoSeq metagenomes, and KO
// detection by other tools (--tools: DIAMOND, fmh-funprofiler, kMermaid, HUMAnN 3 and 4).
// See README.md.

include { MGNIFY_DB; MGNIFY_DB as MGNIFY_MEMBER_DB; MGNIFY_ANNOTATE; MGNIFY_ANNOTATE as MGNIFY_NEAREST_ANNOTATE } from './mgnify_diamond.nf'

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
    # Debian unzip flags this Zenodo archive as overlapping (false-positive zip bomb).
    # Extract aside and move into place only on success: storeDir publishes outputs even
    # when the task fails, so a partial genomes dir would block every later attempt.
    UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE unzip -q -d tmp genomes.zip -x '*.DS_Store'
    mv tmp/genomes_extracted_from_kegg . && rm -r tmp genomes.zip
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

// ---- MGnify90-level truth (--mgnify_indexes with members): every genome protein searched
// against the representatives of the index's clusters (DIAMOND blastp), so detection at the
// 90% level, the nearest cluster beyond it and containment AAI can be checked.

process MGNIFY_INDEX {
    tag "${name}"
    label 'process_high_memory'
    // hard links (same filesystem as the work dir): the query cost run reads them from here
    publishDir "${params.outdir}/mgnify_index", mode: 'link', saveAs: { "${name}/index" }

    input:
    tuple val(name), val(args), path(members, stageAs: 'members/*'), path(pfam, stageAs: 'pfam/*')
    path code, stageAs: 'code/*'

    output:
    tuple val(name), path("index_${name}"), emit: index

    script:
    "${params.kfp} index ${members} index_${name} ${pfam ? "--pfam ${pfam}" : ''} ${args}"

    stub:
    "mkdir index_${name}"
}

process MGNIFY_REPS {
    tag "${name}"
    label 'process_medium'

    input:
    tuple val(name), path(members, stageAs: 'members/*')
    path code, stageAs: 'code/*'

    output:
    tuple val(name), path('reps.faa'), emit: faa

    script:
    "${params.bench} reps --members ${members} --out reps.faa"

    stub:
    "touch reps.faa"
}

process MGNIFY_GENES {
    tag "${name}"
    label 'process_medium'
    publishDir "${params.outdir}/mgnify", mode: 'copy', saveAs: { "${name}_gene_units_reps.parquet" }

    input:
    tuple val(name), path(hits, stageAs: 'hits/*')
    path code, stageAs: 'code/*'

    output:
    tuple val(name), path('gene_units.parquet'), emit: genes

    script:
    "${params.bench} mgnify-genes --hits hits/* --out gene_units.parquet"

    stub:
    "touch gene_units.parquet"
}

process MGNIFY_MEMBERS {
    tag "${name}"
    label 'process_medium'

    input:
    tuple val(name), path(gene_units), path(members, stageAs: 'members/*')
    path code, stageAs: 'code/*'

    output:
    tuple val(name), path('members.faa'), emit: faa
    tuple val(name), path('member_clusters.parquet'), emit: clusters

    script:
    "${params.bench} mgnify-members --members ${members} --gene-units ${gene_units}"

    stub:
    "touch members.faa member_clusters.parquet"
}

process MGNIFY_NEAREST {
    tag "${name}"
    label 'process_high_memory'
    publishDir "${params.outdir}/mgnify", mode: 'copy', saveAs: { "${name}_gene_units.parquet" }

    input:
    tuple val(name), path(hits, stageAs: 'hits/*'), path(gene_units), path(clusters)
    path code, stageAs: 'code/*'

    output:
    tuple val(name), path('gene_units_nearest.parquet'), emit: genes

    script:
    """
    ${params.bench} mgnify-nearest --hits hits/* --gene-units ${gene_units} \\
        --member-clusters ${clusters} --min-cov ${params.mgnify_nearest_min_cov} \\
        --out gene_units_nearest.parquet
    """

    stub:
    "touch gene_units_nearest.parquet"
}

process HOLDOUT_MEMBERS {
    tag "${name} h${fraction}"
    label 'process_medium'

    input:
    tuple val(name), val(fraction), val(args), path(members, stageAs: 'members/*'), path(gene_units)
    path genomes
    path sampled
    path code, stageAs: 'code/*'

    output:
    tuple val(name), val(args), path('held/members.parquet'), emit: members

    script:
    """
    mkdir held
    ${params.bench} holdout-members --members ${members} --gene-units ${gene_units} \
        --genomes-dir ${genomes} --genomes-list ${sampled} --fraction ${fraction} \
        --min-id ${params.holdout_min_id} --out held/members.parquet
    """

    stub:
    "mkdir held && touch held/members.parquet"
}

process UNKNOWN_SCORE {
    tag "seed ${sid} ${name}${arm ? '~' + arm : ''}"
    label 'process_single'

    input:
    tuple val(sid), val(name), val(arm), val(fraction), path(summary), path(gene_units), path(genes)
    path genomes
    path sampled
    path code, stageAs: 'code/*'

    output:
    path 'unknown_score.tsv', emit: score

    script:
    """
    ${params.bench} unknown-score --summary ${summary} --gene-units ${gene_units} --genes ${genes} \
        --genomes-dir ${genomes} --genomes-list ${sampled} --fraction ${fraction} \
        --min-id ${params.holdout_min_id} --sample seed${sid} --index ${name} --arm '${arm}'
    """

    stub:
    "touch unknown_score.tsv"
}

process UNKNOWN_SUMMARY {
    label 'process_single'
    publishDir params.outdir, mode: 'copy'

    input:
    path scores, stageAs: 'unknown*.tsv'

    output:
    path 'unknown_summary.tsv'
    path 'unknown_scores.tsv'

    script:
    "${params.bench} summary ${scores} --keys index arm holdout --out unknown_summary.tsv --scores-out unknown_scores.tsv"

    stub:
    "touch unknown_summary.tsv unknown_scores.tsv"
}

process AAI_SCORE {
    tag "seed ${sid} ${name}${arm ? '~' + arm : ''}${model_name ? '+' + model_name : ''}"
    label 'process_single'

    input:
    // model: '' for the profile's own aai, else an aai_model.json (params.aai_models) to
    // re-estimate it under, scored as arm '<arm>+<model_name>'
    tuple val(sid), val(name), val(arm), path(units), path(gene_units), path(genes), val(model_name), val(model)
    path code, stageAs: 'code/*'

    output:
    path 'aai_score.tsv', emit: score

    script:
    def scored = model_name ? "${arm}+${model_name}" : arm
    """
    ${params.bench} aai-score --profile ${units} --gene-units ${gene_units} --genes ${genes} \\
        --sample seed${sid} --index ${name} --arm '${scored}' --min-id ${params.mgnify_min_id} \\
        --min-cov ${params.mgnify_min_cov} ${model ? "--model ${model}" : ''}
    """

    stub:
    "touch aai_score.tsv"
}

process AAI_CALIBRATE {
    tag "${name}${arm ? '~' + arm : ''}"
    label 'process_single'
    publishDir "${params.outdir}/calibration", mode: 'copy', saveAs: { f -> "${name}${arm ? '~' + arm : ''}${f.endsWith('.json') ? '.json' : '_scores.tsv'}" }

    input:
    // units and genes in the same seed order: the step pairs them by position
    tuple val(name), val(arm), path(units, stageAs: 'units/u*.tsv'), path(genes, stageAs: 'genes/g*.csv'), path(gene_units), path(index)
    path code, stageAs: 'code/*'

    output:
    path 'aai_calibration.json'
    path 'aai_calibration_scores.tsv'

    script:
    """
    ${params.bench} aai-calibrate --profiles ${units} --genes ${genes} --gene-units ${gene_units} \
        --index ${index} --min-id ${params.mgnify_min_id} --min-cov ${params.mgnify_min_cov}
    """

    stub:
    "touch aai_calibration.json aai_calibration_scores.tsv"
}

process AAI_SUMMARY {
    label 'process_single'
    publishDir params.outdir, mode: 'copy'

    input:
    path scores, stageAs: 'aai*.tsv'

    output:
    path 'aai_summary.tsv'
    path 'aai_scores.tsv'

    script:
    "${params.bench} summary ${scores} --keys index arm aai_truth --out aai_summary.tsv --scores-out aai_scores.tsv"

    stub:
    "touch aai_summary.tsv aai_scores.tsv"
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

process GENOME_PROTEINS {
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
    tag "seed ${sid}"
    label 'process_medium'
    publishDir "${params.outdir}/truth", mode: 'copy', saveAs: { f -> f == 'truth.csv' ? "seed${sid}.csv" : "seed${sid}_${f - 'truth_'}" }

    input:
    tuple val(sid), path(fna), path(genes), path(r1), path(r2)
    path kos
    path domains  // [] without Pfam
    path code, stageAs: 'code/*'  // bench.py: only here so -resume reruns on changes

    output:
    tuple val(sid), val('ko'), path('truth.csv'), emit: truth
    tuple val(sid), val('pfam'), path('truth_pfam.csv'), emit: pfam, optional: true
    tuple val(sid), path('truth_genes.csv'), emit: genes

    script:
    """
    ${params.bench} truth --fna ${fna} --genes ${genes} --kos ${kos} --r1 ${r1} --r2 ${r2} \\
        --threads ${task.cpus} ${domains ? "--domains ${domains}" : ''}
    """

    stub:
    "touch truth.csv truth_genes.csv ${domains ? 'truth_pfam.csv' : ''}"
}

process PROFILE {
    tag "seed ${sid} ${name}${arm ? '~' + arm : ''}"
    label 'process_single'
    // one saveAs for both outputs: without the per-file branch they overwrite each other
    publishDir params.outdir, mode: 'copy', saveAs: { f ->
        def id = "seed${sid}_${name}${arm ? '~' + arm : ''}"
        f.endsWith('.parquet') ? "kmers/${id}.parquet" : f == 'units.tsv' ? "units/${id}.tsv" :
            f == 'summary.json' ? "summaries/${id}.json" : "profiles/${id}.tsv"
    }

    input:
    tuple val(sid), path(r1), path(r2), val(label), val(name), path(index), val(by_pfam), val(arm), val(args), path(mask, stageAs: 'mask'), path(decoy, stageAs: 'decoy')
    path code, stageAs: 'code/*'  // query sources: only here so -resume reruns on changes

    output:
    tuple val(sid), val(label), val(name), val(arm), path('profile.tsv'), path('kmers.parquet'), emit: profile
    tuple val(sid), val(name), val(arm), path('units.tsv'), emit: units, optional: true
    tuple val(sid), val(name), val(arm), path('summary.json'), emit: summary, optional: true

    script:
    // by_pfam: units carry Pfam labels (MGnify90 clusters), and the profile is summed per Pfam
    def out = by_pfam ? 'units.tsv' : 'profile.tsv'
    def per_pfam = "${params.bench} pfam-profile --profile units.tsv --unit-pfam ${index}/unit_pfam.parquet"
    def extra = (mask ? ' --mask mask' : '') + (decoy ? ' --extra-index decoy' : '') +
        (params.holdout && by_pfam ? ' --summary summary.json' : '')
    """
    ${params.kfp} query ${index} ${r1} ${r2} --out ${out} --draws ${params.draws} --kmers kmers.parquet \\
        --all-estimators ${args}${extra}
    ${by_pfam ? "${per_pfam} --out profile.tsv" : ''}
    """

    stub:
    "touch profile.tsv kmers.parquet ${by_pfam ? 'units.tsv' : ''}"
}

process SCORE {
    tag "seed ${sid} ${name}${arm ? '~' + arm : ''}"
    label 'process_single'

    input:
    tuple val(sid), val(label), val(name), val(arm), path(profile), path(truth)
    path code, stageAs: 'code/*'  // bench.py: only here so -resume reruns on changes

    output:
    path 'score.tsv', emit: score

    script:
    // DIAMOND's evidence is read pairs, not k-mers: its own range, to trace a curve
    def min_hits = name == 'diamond' ? params.diamond_min_hits : params.min_hits
    """
    ${params.bench} score --truth ${truth} --profile ${profile} --sample seed${sid} --index ${name} \\
        --label ${label} --arm '${arm}' --min-hits ${min_hits.toString().tokenize(',').join(' ')}
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

// ---- Phase-7 arms: read QC (fastp), host spike-in (simulated human reads mixed in at
// constant depth), host handling (hostile upstream, the mask sidecar, a human-proteome decoy
// queried jointly). See README, Ablations.

process FASTP {
    tag "seed ${sid}"
    label 'process_medium'
    container 'quay.io/biocontainers/fastp:1.0.1--heae3180_0'

    input:
    tuple val(sid), path(r1), path(r2)

    output:
    tuple val('fastp'), val(sid), path('fastp_R1.fastq.gz'), path('fastp_R2.fastq.gz'), emit: reads
    path 'fastp.json'

    script:
    """
    fastp -i ${r1} -I ${r2} -o fastp_R1.fastq.gz -O fastp_R2.fastq.gz --trim_poly_g \\
        --thread ${task.cpus} --json fastp.json --html fastp.html ${params.fastp_args}
    """

    stub:
    "touch fastp_R1.fastq.gz fastp_R2.fastq.gz fastp.json"
}

process HOSTILE_DB {
    label 'process_single'
    storeDir "${params.db_dir}/hostile"
    container 'quay.io/biocontainers/hostile:2.0.2--pyhdfd78af_0'

    output:
    path 'cache', emit: db

    script:
    "HOSTILE_CACHE_DIR=\$PWD/cache hostile index fetch --name ${params.hostile_index} --bowtie2"

    stub:
    "mkdir cache"
}

process HOSTILE {
    tag "seed ${sid}"
    label 'process_medium'
    container 'quay.io/biocontainers/hostile:2.0.2--pyhdfd78af_0'

    input:
    tuple val(sid), path(r1), path(r2)
    path db

    output:
    tuple val('hostile'), val(sid), path('clean/*.clean_1.fastq.gz'), path('clean/*.clean_2.fastq.gz'), emit: reads

    script:
    """
    HOSTILE_CACHE_DIR=\$PWD/${db} hostile clean --fastq1 ${r1} --fastq2 ${r2} --aligner bowtie2 \\
        --index ${params.hostile_index} --airplane --threads ${task.cpus} --output clean
    """

    stub:
    "mkdir clean && touch clean/r.clean_1.fastq.gz clean/r.clean_2.fastq.gz"
}

process HOST_GENOME {
    label 'process_single'
    storeDir "${params.db_dir}/host"

    output:
    path 'host.fa.gz', emit: fasta
    path 'abundance.txt', emit: abundance

    script:
    // gzip members concatenate: the genome, then the extra records (rCRS chrM, PhiX)
    """
    curl -fsSL '${params.host_genome_url}' > host.fa.gz
    curl -fsSL '${params.host_extra_url}' | gzip >> host.fa.gz
    ${params.bench} host-abundance --fasta host.fa.gz
    """

    stub:
    "touch host.fa.gz abundance.txt"
}

process HOST_READS {
    label 'process_medium'

    input:
    path fasta
    path abundance
    val n_reads

    output:
    tuple path('host_R1.fastq.gz'), path('host_R2.fastq.gz'), emit: reads

    script:
    // the same read model as the metagenomes; one set of host reads serves every sample
    def iss = params.iss_mode == 'perfect' ? "${params.bench} iss generate --mode perfect" : "${params.venv}/bin/iss generate --model ${params.iss_model}"
    """
    gunzip -c ${fasta} > host.fa
    ${iss} --genomes host.fa --abundance_file ${abundance} --n_reads ${n_reads} --seed 0 \\
        --cpus ${task.cpus} --compress --output host
    rm host.fa
    """

    stub:
    "touch host_R1.fastq.gz host_R2.fastq.gz"
}

process MIX {
    tag "seed ${sid}"
    label 'process_single'

    input:
    tuple val(sid), val(seed), path(r1), path(r2), val(fraction)
    tuple path(host_r1), path(host_r2)
    path code, stageAs: 'code/*'  // bench.py: only here so -resume reruns on changes

    output:
    tuple val(sid), val(seed), path('mixed_R1.fastq.gz'), path('mixed_R2.fastq.gz'), emit: reads

    script:
    """
    ${params.bench} mix --r1 ${r1} --r2 ${r2} --host-r1 ${host_r1} --host-r2 ${host_r2} \\
        --fraction ${fraction} --seed ${seed}
    """

    stub:
    "touch mixed_R1.fastq.gz mixed_R2.fastq.gz"
}

process MASK {
    tag "${name}"
    label 'process_medium'
    publishDir "${params.outdir}/masks", mode: 'copy', pattern: 'mask.json', saveAs: { "${name}.json" }

    input:
    tuple val(label), val(name), path(index)
    path host
    path code, stageAs: 'code/*'  // package sources: only here so -resume rebuilds on changes

    output:
    tuple val(label), val(name), path("mask_${name}"), emit: mask
    path 'mask.json'  // masked hashes and postings (the mask's cost), copied out to publish

    script:
    "${params.kfp} mask ${host} ${index} mask_${name} && cp mask_${name}/mask.json mask.json"

    stub:
    "mkdir mask_${name} && touch mask_${name}/mask.json mask.json"
}

process DECOY_PROTEOME {
    label 'process_single'
    storeDir "${params.db_dir}/decoy"

    input:
    path code, stageAs: 'code/*'

    output:
    path 'decoy_members.parquet', emit: members

    script:
    """
    curl -fsSL '${params.decoy_proteome_url}' > proteome.fa.gz
    ${params.bench} decoy-members --faa proteome.fa.gz
    rm proteome.fa.gz
    """

    stub:
    "touch decoy_members.parquet"
}

process DECOY_INDEX {
    tag "${cfg}"
    label 'process_high_memory'

    input:
    tuple val(cfg), val(args)
    path members
    path code, stageAs: 'code/*'

    output:
    tuple val(cfg), path("decoy_${cfg}"), emit: index

    script:
    // the index config's own arguments, so k, alphabet, hash and tiers match for a joint query
    "${params.kfp} index ${members} decoy_${cfg} ${args}"

    stub:
    "mkdir decoy_${cfg}"
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

// `--tools ''` arrives as true: no tools
def toolList() {
    def tools = params.tools.toString() in ['', 'true', 'false'] ? [] : params.tools.toString().tokenize(',')*.trim()
    def unknown = tools - ['diamond', 'fmh_funprofiler', 'kmermaid', 'humann', 'humann4']
    if (unknown) {
        error "Unknown --tools: ${unknown.join(', ')}"
    }
    return tools
}

// HUMAnN 3.9 (~45 GB) and 4 alpha (~70 GB): ChocoPhlAn full, a UniRef90 DIAMOND database
// (4 distributes only the EC-filtered one), utility mapping, and the MetaPhlAn index each uses.
// Returns the HUMANN_DB and METAPHLAN_DB inputs of the HUMAnN versions in tools.
def humannDbs(tools) {
    def humann = [
        humann:  [params.humann3_container, 'uniref90_diamond', params.metaphlan_index],
        humann4: [params.humann4_container, 'uniref90_ec_filtered_diamond', params.humann4_metaphlan_index],
    ].findAll { tool, _cfg -> tool in tools }
    return [
        humann.collectMany { tool, cfg ->
            [[tool, cfg[0], 'chocophlan', 'full'], [tool, cfg[0], 'uniref', cfg[1]], [tool, cfg[0], 'utility_mapping', 'full']]
        },
        humann.collect { tool, cfg -> [tool, cfg[0], cfg[2]] },
    ]
}

// Databases only: the --tools' databases and the --dbs ones into --db_dir (and the Zenodo
// inputs, which fmh-funprofiler's sketches are, into --data_dir), for later runs with the
// same dirs to reuse
//   nextflow run workflows/fmh-benchmark --dbs_only --tools humann,humann4 --dbs pfam,hostile
workflow DBS {
    def tools = toolList()
    def dbs = params.dbs.toString() in ['', 'true', 'false'] ? [] : params.dbs.toString().tokenize(',')*.trim()
    if (dbs - ['pfam', 'hostile', 'host', 'decoy']) {
        error "Unknown --dbs: ${(dbs - ['pfam', 'hostile', 'host', 'decoy']).join(', ')}"
    }
    if ('pfam' in dbs) {
        PFAM_DB()
    }
    if ('hostile' in dbs) {
        HOSTILE_DB()
    }
    if ('host' in dbs) {
        HOST_GENOME()
    }
    if ('decoy' in dbs) {
        DECOY_PROTEOME(channel.fromPath("${projectDir}/../../python/kmer_functional_profiler/*.py").collect())
    }
    if (tools.any { it in ['diamond', 'fmh_funprofiler', 'kmermaid'] }) {
        FETCH()
    }
    if ('diamond' in tools) {
        DIAMOND_DB(FETCH.out.faa)
    }
    if ('kmermaid' in tools) {
        MEMBERS(FETCH.out.faa, FETCH.out.kos)
        KMERMAID_MODEL(MEMBERS.out.members)
    }
    def (humann_dbs, metaphlan_dbs) = humannDbs(tools)
    HUMANN_DB(channel.fromList(humann_dbs))
    METAPHLAN_DB(channel.fromList(metaphlan_dbs))
}

workflow {
    if (params.dbs_only) {
        DBS()
    } else {
        BENCHMARK()
    }
}

workflow BENCHMARK {
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
    def holdout_fractions = params.holdout.toString().tokenize(',')*.trim()*.toDouble().findAll { it > 0 }
    def ch_domains = channel.value([])
    // MGnify90 indexes: [name, path] (built elsewhere) and/or [members, pfam, args] (built here);
    // members (files or directories) also give the cluster representatives for annotation
    def mgnify = params.mgnify_indexes.collect { cfg ->
        [name: cfg.name, path: cfg.path ?: '', members: cfg.members ?: '', pfam: cfg.pfam ?: '', args: cfg.args ?: '']
    }
    if (mgnify.any { !it.path && !it.members }) {
        error "Each mgnify_indexes entry needs a path or members"
    }
    if ('pfam' in labels || mgnify.any { it.members }) {
        GENOME_PROTEINS(FETCH.out.genomes, bench_py)
    }
    if ('pfam' in labels) {
        def ch_hmm = params.pfam_hmm ? channel.value(file(params.pfam_hmm)) : PFAM_DB().hmm
        PFAM_ANNOTATE(GENOME_PROTEINS.out.faa.flatten(), ch_hmm)
        PFAM_DOMAINS(PFAM_ANNOTATE.out.domtbl.collect(), GENOME_PROTEINS.out.faa, bench_py)
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
    def ch_gene_units = channel.empty()  // name, gene_units
    def ch_rep_units = channel.empty()  // name, gene_units against the representatives
    if (mgnify.any { it.members }) {
        // Entries with the same members (e.g. sparse and dense builds of one subset) share one
        // annotation, run under the first entry's name and handed to the others.
        def sharing = mgnify.findAll { it.members }.groupBy { it.members.toString() }.values()
            .collectEntries { cfgs -> [cfgs[0].name, cfgs*.name] }
        MGNIFY_REPS(
            channel.fromList(mgnify.findAll { it.name in sharing.keySet() }).map { cfg -> [cfg.name, files(cfg.members)] },
            bench_py,
        )
        MGNIFY_DB(MGNIFY_REPS.out.faa)
        MGNIFY_ANNOTATE(MGNIFY_DB.out.db.combine(GENOME_PROTEINS.out.faa.flatten()), params.mgnify_diamond_args)
        MGNIFY_GENES(MGNIFY_ANNOTATE.out.hits.groupTuple(), bench_py)
        // nearest-member identity, what aai estimates (plan, phase 7, step 34): the genes
        // against the members of the clusters they hit
        MGNIFY_MEMBERS(
            MGNIFY_GENES.out.genes.join(
                channel.fromList(mgnify.findAll { it.name in sharing.keySet() }).map { cfg -> [cfg.name, files(cfg.members)] }
            ),
            bench_py,
        )
        MGNIFY_MEMBER_DB(MGNIFY_MEMBERS.out.faa)
        MGNIFY_NEAREST_ANNOTATE(
            MGNIFY_MEMBER_DB.out.db.combine(GENOME_PROTEINS.out.faa.flatten()),
            params.mgnify_nearest_diamond_args,
        )
        MGNIFY_NEAREST(
            MGNIFY_NEAREST_ANNOTATE.out.hits.groupTuple().join(MGNIFY_GENES.out.genes).join(MGNIFY_MEMBERS.out.clusters),
            bench_py,
        )
        ch_gene_units = MGNIFY_NEAREST.out.genes.flatMap { name, genes -> sharing[name].collect { n -> [n, genes] } }
        ch_rep_units = MGNIFY_GENES.out.genes.flatMap { name, genes -> sharing[name].collect { n -> [n, genes] } }
    }
    // Unknown-fraction hold-out ladder (plan, phase 7, step 36): each entry built here is also
    // built without the units the held-out genomes hit, as <name>_h<pct>; holdout_of maps every
    // index in the ladder to [entry, fraction] (the entry itself is fraction 0).
    def holdout_of = [:]
    SAMPLE(channel.of(1..params.replicates), FETCH.out.genomes)
    // genomes the samples drew: the ladder holds out among these only (step 39: from the
    // whole pool it took most shared units while holding out few sampled genomes)
    def ch_sampled = SAMPLE.out[1].collectFile(name: 'sampled_genomes.txt', sort: true).first()
    if (params.holdout) {
        def built = mgnify.findAll { it.members && !it.path }
        if (!built || built.any { !it.pfam }) {
            error "--holdout needs mgnify_indexes entries built here, with members, pfam and args"
        }
        built.each { cfg ->
            holdout_of[cfg.name] = [cfg.name, 0d]
            holdout_fractions.each { h -> holdout_of["${cfg.name}_h${Math.round(h * 100)}".toString()] = [cfg.name, h] }
        }
        HOLDOUT_MEMBERS(
            channel.fromList(built).flatMap { cfg ->
                holdout_fractions.collect { h ->
                    [cfg.name, "${cfg.name}_h${Math.round(h * 100)}".toString(), h, cfg.args, files(cfg.members)]
                }
            }
                .combine(ch_rep_units, by: 0)
                .map { _base, hname, h, args, members, gene_units -> [hname, h, args, members, gene_units] },
            FETCH.out.genomes,
            ch_sampled,
            bench_py,
        )
    }
    def ch_mgnify = channel.fromList(mgnify.findAll { it.path }).map { cfg -> [cfg.name, file(cfg.path, checkIfExists: true)] }
    if (mgnify.any { !it.path }) {
        MGNIFY_INDEX(
            channel.fromList(mgnify.findAll { !it.path }).map { cfg ->
                [cfg.name, cfg.args, files(cfg.members), cfg.pfam ? file(cfg.pfam, checkIfExists: true) : []]
            }.mix(
                // the entry's Pfam table as is: the build joins it on the members kept
                params.holdout
                    ? HOLDOUT_MEMBERS.out.members.map { hname, args, members ->
                        [hname, args, members, file(mgnify.find { it.name == holdout_of[hname][0] }.pfam)]
                    }
                    : channel.empty()
            ),
            ch_code,
        )
        ch_mgnify = ch_mgnify.mix(MGNIFY_INDEX.out.index)
    }
    // MGnify indexes are profiled whatever the labels: their unit profiles carry the AAI truth
    ch_indexes = ch_indexes.mix(ch_mgnify.map { name, index -> ['pfam', name, index, true] })

    // Query arms (phase 7): which reads (raw, fastp, hostile), extra query options, and
    // whether the host mask and the human-proteome decoy are used. Default: one plain arm.
    def arms = params.query_arms.collect { a ->
        [reads: a.reads ?: 'raw', name: a.name ?: '', args: a.args ?: '', mask: a.mask ?: false, decoy: a.decoy ?: false]
    }
    def kinds = arms*.reads.unique()
    if (kinds - ['raw', 'fastp', 'hostile']) {
        error "Unknown query_arms reads: ${(kinds - ['raw', 'fastp', 'hostile']).join(', ')}"
    }
    def fractions = params.host_fractions.toString().tokenize(',')*.trim()*.toDouble()
    def host_fractions = fractions.findAll { it > 0 }

    SIMULATE(SAMPLE.out.sample.map { seed, fna, _genes -> [seed, fna] })
    // sid: the sample's id, the seed, or seed + host percentage (1h90) for spike-in samples
    def ch_reads = SIMULATE.out.reads.filter { 0d in fractions }.map { seed, r1, r2 -> [seed.toString(), seed, r1, r2] }
    if (host_fractions || arms.any { it.mask }) {
        HOST_GENOME()
    }
    if (host_fractions) {
        // host pairs: the largest share of the metagenomes' pairs, with 2% to spare
        def n_host = Math.ceil(params.n_reads * host_fractions.max() * 1.02).toLong() + 100
        HOST_READS(HOST_GENOME.out.fasta, HOST_GENOME.out.abundance, n_host)
        MIX(
            SIMULATE.out.reads.combine(channel.fromList(host_fractions)).map { seed, r1, r2, f ->
                ["${seed}h${Math.round(f * 100)}".toString(), seed, r1, r2, f]
            },
            HOST_READS.out.reads.first(),
            bench_py,
        )
        ch_reads = ch_reads.mix(MIX.out.reads)
    }
    TRUTH(
        ch_reads.map { sid, seed, r1, r2 -> [seed, sid, r1, r2] }
            .combine(SAMPLE.out.sample, by: 0)
            .map { _seed, sid, r1, r2, fna, genes -> [sid, fna, genes, r1, r2] },
        FETCH.out.kos,
        ch_domains,
        bench_py,
    )
    // sid, label, truth; MGnify profiles (label pfam) are scored on Pfam only with that label
    ch_truth = TRUTH.out.truth.mix(TRUTH.out.pfam).filter { it[1] in labels }

    def ch_sid_reads = ch_reads.map { sid, _seed, r1, r2 -> [sid, r1, r2] }
    def ch_variants = ch_sid_reads.map { sid, r1, r2 -> ['raw', sid, r1, r2] }
    if ('fastp' in kinds) {
        FASTP(ch_sid_reads)
        ch_variants = ch_variants.mix(FASTP.out.reads)
    }
    if ('hostile' in kinds) {
        HOSTILE(ch_sid_reads, HOSTILE_DB().db)
        ch_variants = ch_variants.mix(HOSTILE.out.reads)
    }
    // Masks and decoys exist for the kfp-hashed indexes built here (not fmh_compat or MGnify)
    def ch_masks = channel.empty()
    if (arms.any { it.mask }) {
        MASK(INDEX.out.index, HOST_GENOME.out.fasta, ch_code)
        ch_masks = MASK.out.mask
    }
    def ch_decoys = channel.empty()
    if (arms.any { it.decoy }) {
        DECOY_INDEX(
            channel.fromList(params.indexes).map { cfg -> [cfg.name, cfg.args] },
            DECOY_PROTEOME(ch_code).members,
            ch_code,
        )
        ch_decoys = DECOY_INDEX.out.index
    }
    // label, name, index, by_pfam, mask (or null), decoy (or null)
    def ch_extended = ch_indexes
        .join(ch_masks, by: [0, 1], remainder: true)
        .map { label, name, index, by_pfam, mask -> [label == 'ko' ? name : name - "${label}_", label, name, index, by_pfam, mask] }
        .combine(ch_decoys.ifEmpty(['', null]).toList().map { it.collectEntries() })
        .map { cfg, label, name, index, by_pfam, mask, decoys -> [label, name, index, by_pfam, mask, decoys[cfg]] }
        .filter { it[2] != null }  // join's remainder of an index without a mask
    def ch_runs = ch_variants
        .combine(channel.fromList(arms).map { a -> [a.reads, a] }, by: 0)
        .combine(ch_extended)
        .filter { _kind, _sid, _r1, _r2, arm, _label, _name, _index, _by_pfam, mask, decoy ->
            (!arm.mask || mask != null) && (!arm.decoy || decoy != null)
        }
        .map { _kind, sid, r1, r2, arm, label, name, index, by_pfam, mask, decoy ->
            [sid, r1, r2, label, name, index, by_pfam, arm.name, arm.args, arm.mask ? mask : [], arm.decoy ? decoy : []]
        }

    PROFILE(ch_runs, ch_code)
    // sid, label, name, arm, profile, kmers, truth
    ch_scored = PROFILE.out.profile.combine(ch_truth, by: [0, 1])
    DETECTED(
        ch_scored
            .filter { it[3] == '' }  // diagnostics for the plain arm only
            .combine(ch_reads.map { sid, seed, _r1, _r2 -> [seed, sid] }.combine(SAMPLE.out.sample, by: 0)
                .map { _seed, sid, fna, _genes -> [sid, fna] }, by: 0)
            .map { sid, label, name, _arm, profile, kmers, truth, fna -> [label, name, sid, profile, kmers, truth, fna] }
            .combine(ch_indexes.filter { !it[3] }.map { it.take(3) }, by: [0, 1]),
        bench_py,
    )
    DETECTED.out.detected.collectFile(name: 'detected.tsv', keepHeader: true, storeDir: params.outdir)

    def tools = toolList()
    def ch_raw = channel.empty()
    def ch_tool_reads = SIMULATE.out.reads.filter { 0d in fractions }.map { seed, r1, r2 -> [seed.toString(), r1, r2] }
    if ('diamond' in tools) {
        DIAMOND(ch_tool_reads, DIAMOND_DB(FETCH.out.faa).db)
        ch_raw = ch_raw.mix(DIAMOND.out.raw)
    }
    if ('fmh_funprofiler' in tools) {
        FMH_FUNPROFILER(ch_tool_reads, FETCH.out.sketches)
        ch_raw = ch_raw.mix(FMH_FUNPROFILER.out.raw)
    }
    if ('kmermaid' in tools) {
        KMERMAID(ch_tool_reads, KMERMAID_MODEL(MEMBERS.out.members).model)
        ch_raw = ch_raw.mix(KMERMAID.out.raw)
    }
    def dbs = ['chocophlan', 'uniref', 'utility_mapping']
    def (humann_dbs, metaphlan_dbs) = humannDbs(tools)
    HUMANN_DB(channel.fromList(humann_dbs))
    METAPHLAN_DB(channel.fromList(metaphlan_dbs))
    // tool -> [chocophlan, uniref, utility_mapping, metaphlan], as HUMANN and HUMANN4 take them
    def ch_humann_db = HUMANN_DB.out.db
        .groupTuple(size: dbs.size())
        .join(METAPHLAN_DB.out.db)
        .map { tool, names, paths, mpa -> [tool] + dbs.collect { db -> paths[names.indexOf(db)] } + [mpa] }
    if ('humann' in tools) {
        HUMANN(ch_tool_reads, ch_humann_db.filter { it[0] == 'humann' }.map { it.drop(1) }.first())
        ch_raw = ch_raw.mix(HUMANN.out.raw)
    }
    if ('humann4' in tools) {
        HUMANN4(ch_tool_reads, ch_humann_db.filter { it[0] == 'humann4' }.map { it.drop(1) }.first())
        ch_raw = ch_raw.mix(HUMANN4.out.raw)
    }
    TOOL_PROFILE(ch_raw, FETCH.out.kos, MEMBERS.out.members, bench_py)

    SCORE(
        ch_scored
            .map { sid, label, name, arm, profile, _kmers, truth -> [sid, label, name, arm, profile, truth] }
            .mix(TOOL_PROFILE.out.profile.map { sid, tool, profile -> [sid, 'ko', tool, '', profile] }
                .combine(TRUTH.out.truth.map { sid, label, truth -> [sid, label, truth] }, by: [0, 1])
                .map { sid, label, tool, arm, profile, truth -> [sid, label, tool, arm, profile, truth] }),
        bench_py,
    )
    SUMMARY(SCORE.out.score.collect())
    // unit-level detection and AAI of the MGnify indexes against the DIAMOND truth
    def ch_aai = PROFILE.out.units
        .map { sid, name, arm, units -> [name, sid, arm, units] }
        .combine(ch_gene_units, by: 0)  // name, sid, arm, units, gene_units
        .map { name, sid, arm, units, gene_units -> [sid, name, arm, units, gene_units] }
        .combine(TRUTH.out.genes, by: 0)  // sid, name, arm, units, gene_units, genes
    def ch_models = channel.of(['', '']).mix(
        channel.fromList(params.aai_models).map { m -> [m.name, file(m.path, checkIfExists: true).toString()] }
    )
    AAI_SCORE(ch_aai.combine(ch_models), bench_py)
    AAI_SUMMARY(AAI_SCORE.out.score.collect())
    if (params.holdout) {
        // sid, name, arm, fraction, summary, gene_units (the entry's), truth genes
        UNKNOWN_SCORE(
            PROFILE.out.summary
                .filter { _sid, name, _arm, _summary -> holdout_of.containsKey(name) }
                .map { sid, name, arm, summary -> [holdout_of[name][0], sid, name, arm, holdout_of[name][1], summary] }
                .combine(ch_rep_units, by: 0)
                .map { _base, sid, name, arm, h, summary, gene_units -> [sid, name, arm, h, summary, gene_units] }
                .combine(TRUTH.out.genes, by: 0),
            FETCH.out.genomes,
            ch_sampled,
            bench_py,
        )
        UNKNOWN_SUMMARY(UNKNOWN_SCORE.out.score.collect())
    }
    // aai -> identity map per index and arm, fitted on half the clusters of every seed
    AAI_CALIBRATE(
        ch_aai
            .map { sid, name, arm, units, gene_units, genes -> [[name, arm], sid, units, genes, gene_units] }
            .groupTuple()
            .map { key, sids, units, genes, gene_units ->
                def order = (0..<sids.size()).toList().sort { sids[it] }
                [key[0], key[1], order.collect { units[it] }, order.collect { genes[it] }, gene_units[0]]
            }
            .combine(ch_mgnify, by: 0),  // name, arm, units, genes, gene_units, index
        bench_py,
    )

    // run.json: what produced the results in outdir (reads, indexes, code, status). params and
    // workflow are read here: inside the handler, names resolve against the workflow metadata.
    def run_params = params
    def wf = workflow
    def reads = params.iss_mode == 'perfect' ? 'error-free (iss perfect)' : "iss ${params.iss_model} (${params.iss_mode})"
    def description = "${params.replicates} metagenomes x ${params.n_genomes} KEGG genomes, " +
        "${params.n_reads} reads, ${reads}; labels ${params.labels}; ${params.indexes.size()} index configs" +
        " + fmh_compat + ${params.mgnify_indexes.size()} MGnify; ${params.query_arms.size()} query arms;" +
        " host fractions ${params.host_fractions}; draws ${params.draws}"
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
