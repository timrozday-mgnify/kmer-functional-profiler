// Read recruitment and crude pileups downstream of a kfp profile (plan, phase 13; the
// nucleotide follow-up of phase 14, step 4). Per sample: kfp query -> per-sample reference
// (kfp-pileup reference) -> read routing by k-mers to blocks of the reference (kfp-pileup
// route) -> DIAMOND blastx per block (amino-acid index, frameshift-aware) or bwa-mem2
// mem (nucleotide index), chosen by the index's molecule -> read posteriors, pileups,
// per-unit identity, read assignment and, for DIAMOND, the nucleotide pileup and consensus
// (kfp-pileup call) -> the hand-off BAM (reads re-mapped to the consensus, or the bwa-mem2
// alignments, kept where step 8 assigns them) -> variants, consequences, haplotypes and
// microdiversity by existing tools (plan, phase 13, steps 6b-6c). See README.md.

process INDEX {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'

    input:
    path targets
    val opts  // options as inputs, not params in the script: Nextflow does not hash params

    output:
    path 'index', emit: index

    script:
    "${params.kfp} index ${targets} index --targets ${opts}"

    stub:
    "mkdir index && echo {} > index/meta.json"
}

// kfp's query code, a QUERY input: without it -resume kept a profile from before a presence
// change (pileup-cost's cached QUERY predates presence-fdr)
def query_code() {
    ['query', 'partition'].collect {
        file("${moduleDir}/../../python/kmer_functional_profiler/${it}.py", checkIfExists: true)
    }
}

process QUERY {
    tag "${sample}"
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/${it}" }

    input:
    tuple val(sample), path(r1), path(r2)
    path index
    path meta, stageAs: 'index.meta.json'  // see PILEUP
    val opts  // e.g. '--stats stats.json' (query_stats), --aai, --min-identity and query_args
    path code, stageAs: 'code/*'  // unused: staged so that -resume reruns the query when its code changes

    output:
    tuple val(sample), path('profile.tsv'), emit: profile
    path 'stats.json', optional: true  // with --query_stats: stage times, hits per pair

    script:
    "${params.kfp} query ${index} ${r1} ${r2} --out profile.tsv ${opts}"

    stub:
    "touch profile.tsv"
}

process REFERENCE {
    tag "${sample}"
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/${it}" }

    input:
    tuple val(sample), path(profile)
    path index
    path meta, stageAs: 'index.meta.json'  // see PILEUP
    path store
    path withhold  // units left out of the reference (spike-benchmark), or NO_WITHHOLD
    val opts
    path code  // unused: staged so that -resume rebuilds the reference when its code changes

    output:
    tuple val(sample), path('ref'), env('MOLECULE'), emit: ref

    script:
    def stores = store.collect { "--store ${it}" }.join(' ')
    def held = withhold.name == 'NO_WITHHOLD' ? '' : "--withhold ${withhold}"
    """
    ${params.kfp_pileup} reference ${profile} ${index} ref ${opts} ${stores} ${held}
    MOLECULE=\$(cat ref/molecule.txt)
    """

    stub:
    "mkdir ref && echo protein > ref/molecule.txt && MOLECULE=protein"
}

process ROUTE {
    tag "${sample}"
    label 'process_single'
    publishDir params.outdir, mode: 'copy', pattern: 'route*.{json,parquet}', saveAs: { "${sample}/${it}" }

    input:
    tuple val(sample), path(ref), path(r1), path(r2)
    path meta, stageAs: 'index/meta.json'  // the index's k and alphabet
    val opts
    path code  // unused: staged so that -resume reruns routing when its code changes

    output:
    tuple val(sample), path('route'), path(r1), path(r2), emit: route
    path 'route.json'  // pairs routed and dropped, blocks
    path 'routes.parquet'  // every (pair, block) hit and its distinct k-mers

    script:
    """
    ${params.kfp_pileup} route ${ref} index ${r1} ${r2} route ${opts}
    cp route/route.json route/routes.parquet .
    """

    stub:
    "mkdir route && echo 1 > route/letters.txt && echo {} > route.json && touch routes.parquet"
}

process DIAMOND {
    tag "${sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/diamond:2.2.8--he361c42_0'

    input:
    tuple val(sample), path(ref), path(r1), path(r2)
    val opts

    output:
    tuple val(sample), val('diamond'), path('r1.tsv'), path('r2.tsv'), emit: alignments

    script:
    // every hit within --top % of the read's best: the posteriors need the runners-up.
    // Routed (ROUTE's letters.txt): each block's pairs against that block only, e-values
    // on the whole reference's size; a pair in two blocks is sorted back to one run of
    // lines per read (CALL reads a read's alignments as one run)
    """
    align() {  # fasta, mate 1, mate 2, extra options
        diamond makedb --in \$1 --db ref --threads ${task.cpus}
        diamond blastx --db ref.dmnd --query \$2 --threads ${task.cpus} ${opts} \$4 >> r1.tsv
        diamond blastx --db ref.dmnd --query \$3 --threads ${task.cpus} ${opts} \$4 >> r2.tsv
    }
    : > r1.tsv
    : > r2.tsv
    if [ -e ${ref}/letters.txt ]; then
        for q1 in ${ref}/b*_1.fa.gz; do
            [ -e "\$q1" ] || continue
            b=\${q1%_1.fa.gz}
            align \$b.fa \$q1 \${b}_2.fa.gz "--dbsize \$(cat ${ref}/letters.txt)"
        done
        for m in 1 2; do
            LC_ALL=C sort -s -t "\$(printf '\\t')" -k1,1 -o r\$m.tsv r\$m.tsv
        done
    else
        align ${ref}/ref.fa ${r1} ${r2} ''
    fi
    """

    stub:
    "touch r1.tsv r2.tsv"
}

process BWAMEM2 {
    tag "${sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he70b90d_8'

    input:
    tuple val(sample), path(ref), path(r1), path(r2)
    val opts

    output:
    tuple val(sample), val('bwa'), path('aln.sam'), emit: alignments

    script:
    // -a: every alignment, secondaries included, for the posteriors
    """
    bwa-mem2 index -p ref ${ref}/ref.fa
    bwa-mem2 mem -a -t ${task.cpus} ${opts} ref ${r1} ${r2} > aln.sam
    """

    stub:
    "touch aln.sam"
}

process CALL {
    tag "${sample}"
    label 'process_medium'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/${it}" }

    input:
    tuple val(sample), val(aligner), path(alignments), path(ref), path(profile), path(r1), path(r2)
    val opts
    path code, stageAs: 'code/*'  // unused: staged so that -resume reruns CALL when its code changes

    output:
    tuple val(sample), path('pileup'), emit: pileup
    tuple val(sample), val(aligner), path('pileup'), emit: called
    path 'call_stats.json'

    script:
    // DIAMOND: the reads too, for the nucleotide pileup and consensus (step 6b); per-stage
    // time and memory in call_stats.json (CALL's memory, phase 13)
    def reads = aligner == 'diamond' ? "--reads ${r1} --reads ${r2}" : ''
    """
    ${params.kfp_pileup} call ${ref} pileup ${alignments} --aligner ${aligner} \\
        --profile ${profile} ${reads} --stats call_stats.json --threads ${task.cpus} ${opts}
    """

    stub:
    "mkdir pileup; touch call_stats.json"
}

process MAP_CONSENSUS {
    tag "${sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he70b90d_8'

    input:
    tuple val(sample), path(pileup)
    val opts

    output:
    tuple val(sample), path('consensus.sam'), path('consensus.fna'), emit: sam

    script:
    // the recruited pairs re-mapped to the protein arm's consensus coding sequences (6b)
    """
    cp ${pileup}/consensus.fna .
    bwa-mem2 index consensus.fna
    bwa-mem2 mem -t ${task.cpus} ${opts} consensus.fna ${pileup}/recruited_1.fq.gz \\
        ${pileup}/recruited_2.fq.gz > consensus.sam
    """

    stub:
    "touch consensus.sam consensus.fna"
}

process HANDOFF {
    tag "${sample}"
    label 'process_single'

    input:
    tuple val(sample), path(sam), path(fasta), path(pileup), path(ref)
    val opts
    path code  // unused: staged so that -resume reruns the hand-off when its code changes

    output:
    tuple val(sample), path('handoff.sam'), path(fasta), emit: sam

    script:
    "${params.kfp_pileup} handoff ${pileup} ${ref} ${sam} --out handoff.sam ${opts}"

    stub:
    "touch handoff.sam"
}

process BAM {
    tag "${sample}"
    label 'process_single'
    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/variants/${it}" }

    input:
    tuple val(sample), path(sam), path(fasta)

    output:
    tuple val(sample), path('handoff.bam'), path('handoff.bam.bai'), path(fasta), path("${fasta}.fai"), emit: bam

    script:
    """
    samtools sort -o handoff.bam ${sam}
    samtools index handoff.bam
    samtools faidx ${fasta}
    """

    stub:
    "touch handoff.bam handoff.bam.bai ${fasta}.fai"
}

process FREEBAYES {
    tag "${sample}"
    label 'process_single'
    container 'quay.io/biocontainers/freebayes:1.3.10--h3752d28_1'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/variants/${it}" }

    input:
    tuple val(sample), path(bam), path(bai), path(fasta), path(fai)
    val opts

    output:
    tuple val(sample), path('freebayes.vcf.gz'), path('freebayes.vcf.gz.tbi'), emit: vcf

    script:
    // allele frequencies without a ploidy; SNVs, indels and MNPs in one haplotype-based pass
    """
    freebayes -f ${fasta} --pooled-continuous ${opts} ${bam} | bgzip > freebayes.vcf.gz
    tabix -p vcf freebayes.vcf.gz
    """

    stub:
    "touch freebayes.vcf.gz freebayes.vcf.gz.tbi"
}

process LOFREQ {
    tag "${sample}"
    label 'process_single'
    container 'quay.io/biocontainers/lofreq:2.1.5--py311h6385618_16'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/variants/${it}" }

    input:
    tuple val(sample), path(bam), path(bai), path(fasta), path(fai)
    val opts

    output:
    tuple val(sample), path('lofreq.vcf.gz'), emit: vcf

    script:
    // low-frequency SNVs and indels: indel qualities (Dindel), then its Poisson-binomial calls
    """
    lofreq indelqual --dindel -f ${fasta} -o indelqual.bam ${bam}
    samtools index indelqual.bam
    lofreq call --call-indels -f ${fasta} -o lofreq.vcf ${opts} indelqual.bam
    bgzip lofreq.vcf
    """

    stub:
    "touch lofreq.vcf.gz"
}

process CONSEQUENCES {
    tag "${sample}"
    label 'process_single'
    container 'quay.io/biocontainers/bcftools:1.24--h118bc1c_4'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/variants/${it}" }

    input:
    tuple val(sample), path(vcf), path(tbi), path(fasta), path(fai), path(gff)

    output:
    tuple val(sample), path('consequences.vcf.gz'), path('consequences.tsv'), emit: csq

    script:
    // each consensus sequence is one CDS in the protein's frame; local consequences per
    // record (pooled calls have no meaningful genotypes); BCSQ's amino-acid positions are
    // the reference protein's columns (1-based), as the consensus keeps its frame
    """
    bcftools csq -f ${fasta} -g ${gff} --samples - --local-csq -Oz -o consequences.vcf.gz ${vcf}
    printf 'seq_id\\tpos\\tref\\talt\\tqual\\tdepth\\talt_depth\\tconsequence\\n' > consequences.tsv
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/DP\\t%INFO/AO\\t%INFO/BCSQ\\n' \\
        consequences.vcf.gz >> consequences.tsv
    """

    stub:
    "touch consequences.vcf.gz consequences.tsv"
}

process FLORIA {
    tag "${sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/floria:0.0.2--ha6fb395_0'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/variants/${it}" }

    input:
    tuple val(sample), path(bam), path(bai), path(fasta), path(fai), path(vcf), path(tbi)
    val opts

    output:
    tuple val(sample), path('floria'), emit: haplotypes

    script:
    // Floria 0.0.2 panics on a VCF without records: nothing to phase
    """
    if zcat ${vcf} | grep -qv '^#'; then
        floria -b ${bam} -v ${vcf} -r ${fasta} -o floria -t ${task.cpus} ${opts}
    else
        mkdir floria
    fi
    """

    stub:
    "mkdir floria"
}

process CLIQUESNV {
    tag "${sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/cliquesnv:2.0.3--hdfd78af_0'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/variants/${it}" }

    input:
    tuple val(sample), path(sam), path(fasta)
    val min_reads
    val opts

    output:
    tuple val(sample), path('cliquesnv'), emit: haplotypes

    script:
    // CliqueSNV phases one reference per input: the hand-off split per sequence (records
    // sorted by reference, each file with its own @SQ line), run where it has at least
    // min_reads records; haplotypes.tsv gathers every sequence's haplotypes and frequencies
    """
    mkdir -p cliquesnv split
    grep -v '^@' ${sam} | sort -t "\$(printf '\\t')" -k3,3 -s | awk -F '\\t' -v OFS='\\t' -v min=${min_reads} '
        function flush() {
            if (n < min) return
            f = "split/" c ".sam"
            print sq[c] > f
            for (i = 1; i <= n; i++) print r[i] > f
            close(f)
        }
        FNR == NR { if (\$1 == "@SQ") { split(\$2, a, ":"); sq[a[2]] = \$0 }; next }
        \$3 != c { flush(); c = \$3; n = 0 }
        # a mate on another sequence: mate unmapped here (each file has one @SQ)
        \$7 != "=" && \$7 != c { \$7 = "*"; \$8 = 0; \$9 = 0; \$2 = or(\$2, 8) }
        { r[++n] = \$0 }
        END { flush() }' ${sam} -
    printf 'seq_id\\thaplotype\\tfrequency\\n' > cliquesnv/haplotypes.tsv
    for f in split/*.sam; do
        [ -e "\$f" ] || continue
        c=\$(basename "\$f" .sam)
        cliquesnv -m snv-illumina -in "\$f" -outDir "cliquesnv/\$c" -threads ${task.cpus} ${opts} \\
            > "cliquesnv/\$c.log" 2>&1 || echo "CliqueSNV failed on \$c: see cliquesnv/\$c.log" >&2
        # no FASTA when it finds no haplotype (too low coverage)
        [ -e "cliquesnv/\$c/\$c.fasta" ] || continue
        grep '>' "cliquesnv/\$c/\$c.fasta" | sed 's/^>//; s/_fr_/\\t/' \\
            | awk -v c="\$c" '{ print c "\\t" \$0 }' >> cliquesnv/haplotypes.tsv
    done
    """

    stub:
    "mkdir cliquesnv"
}

process INSTRAIN {
    tag "${sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/instrain:1.10.0--pyhdfd78af_0'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/variants/${it}" }

    input:
    tuple val(sample), path(bam), path(bai), path(fasta), path(fai), path(genes)
    val opts

    output:
    tuple val(sample), path('instrain'), emit: profile

    script:
    // per-gene microdiversity (nucleotide diversity, linkage) on the protein arm's consensus
    def g = genes.name == 'NO_GENES' ? '' : "-g ${genes}"
    "inStrain profile ${bam} ${fasta} -o instrain -p ${task.cpus} ${g} --skip_plot_generation ${opts}"

    stub:
    "mkdir instrain"
}

// A boolean param: Nextflow 26 keeps `--floria false` from the command line as the string
// "false", which Groovy reads as true; a config's `false` stays a Boolean.
def flag(value) {
    value instanceof String ? value.toBoolean() : value as boolean
}

// Reusable by other pipelines (workflows/read-xref, spike-benchmark): reads (sample, r1, r2),
// the index, the store's members tables (a list, possibly empty) and a file of units to
// withhold from the reference (NO_WITHHOLD for none).
workflow PILEUP {
    take:
    reads
    index
    store
    withhold

    main:
    // The index's meta.json as an input of its own: a directory's cache key is its path, size
    // and mtime, which an in-place rebuild keeps, so -resume reused a profile of the old
    // index (pileup-cost REFERENCE: unit ids not in the index). A rebuild rewrites meta.json.
    def meta = index.map { it.resolve('meta.json') }
    // aai/ani optional (targets need neither): --aai only when asked for, or with theta
    def query_opts = [
        flag(params.query_stats) ? '--stats stats.json' : '',
        flag(params.aai) ? '--aai' : '',
        params.min_identity ? "--min-identity ${params.min_identity}" : '',
        params.query_args,
    ].findAll { it }.join(' ')
    def profiles = QUERY(reads, index, meta, query_opts, query_code()).profile
    def refs = REFERENCE(profiles, index, meta, store, withhold,
        "--max-fdr ${params.max_fdr} --max-targets ${params.max_targets} --neighbours ${params.neighbours} --members ${params.members} --route-members ${params.route_members}",
        file("${moduleDir}/../../python/kfp_pileup/reference.py", checkIfExists: true)).ref
    def to_align = refs.join(reads).branch { sample, ref, molecule, r1, r2 ->
        nucleotide: molecule == 'nucleotide'
            return tuple(sample, ref, r1, r2)
        protein: true
            return tuple(sample, ref, r1, r2)
    }
    def frameshift = params.frameshift ? "-F ${params.frameshift} " : ''
    def diamond = "--outfmt 6 ${params.diamond_fields} --top ${params.diamond_top} ${frameshift}${params.diamond_args}"
    // routing (phase 13, step 3b's "Afterwards"): each pair against the reference blocks
    // its k-mers hit, pairs with no hit dropped; amino-acid references only
    def to_diamond = flag(params.route)
        ? ROUTE(to_align.protein, meta, "--block-letters ${params.route_block_letters} --k ${params.route_k} --min-kmers ${params.route_min_kmers}",
            file("${moduleDir}/../../python/kfp_pileup/route.py", checkIfExists: true)).route
        : to_align.protein
    def alignments = DIAMOND(to_diamond, diamond).alignments
        .map { sample, aligner, r1, r2 -> tuple(sample, aligner, [r1, r2]) }
        .mix(BWAMEM2(to_align.nucleotide, params.bwa_args).alignments)
    def ref_dirs = refs.map { sample, ref, molecule -> tuple(sample, ref) }
    def call_code = ['pileup', 'refine', 'spill', 'nucleotide'].collect {
        file("${moduleDir}/../../python/kfp_pileup/${it}.py", checkIfExists: true)
    }
    def called = CALL(alignments.join(ref_dirs).join(profiles).join(reads),
        "--min-depth ${params.min_depth} --min-freq ${params.min_freq} --novel ${params.novel} " +
        "--top ${params.call_top} --max-bits-below ${params.call_max_bits_below} " +
        "--unit-bits-below ${params.call_unit_bits_below} --min-identity ${params.call_min_identity}",
        call_code)
    def pileups = called.pileup

    // steps 6b-6c: the hand-off BAM, then existing callers on it
    def protein = called.called.filter { it[1] == 'diamond' }.map { s, a, p -> tuple(s, p) }
    def nucleotide_sams = alignments.filter { it[1] == 'bwa' }.join(ref_dirs)
        .map { s, a, sam, ref -> tuple(s, sam, ref.resolve('ref.fa')) }
    def to_handoff = MAP_CONSENSUS(protein, params.bwa_args).sam.mix(nucleotide_sams)
        .join(pileups).join(ref_dirs)
    def handoff = HANDOFF(to_handoff, "--min-posterior ${params.min_assign}",
        file("${moduleDir}/../../python/kfp_pileup/handoff.py", checkIfExists: true)).sam
    def bams = BAM(handoff).bam
    if (flag(params.cliquesnv)) {
        CLIQUESNV(handoff, params.cliquesnv_min_reads, params.cliquesnv_args)
    }
    def vcfs = FREEBAYES(bams, params.freebayes_args).vcf
    if (flag(params.lofreq)) {
        LOFREQ(bams, params.lofreq_args)
    }
    // consequences on the protein arm's consensus (a CDS by construction)
    CONSEQUENCES(vcfs.join(bams).join(protein)
        .map { s, vcf, tbi, bam, bai, fasta, fai, p -> tuple(s, vcf, tbi, fasta, fai, p.resolve('consensus.gff')) })
    if (flag(params.floria)) {
        FLORIA(bams.join(vcfs), params.floria_args)
    }
    if (flag(params.instrain)) {
        def no_genes = file("${projectDir}/assets/NO_GENES")
        def genes = called.called.map { s, a, p ->
            tuple(s, a == 'diamond' ? p.resolve('consensus.genes.fna') : no_genes) }
        INSTRAIN(bams.join(genes), params.instrain_args)
    }

    emit:
    profile = profiles
    ref = ref_dirs
    pileup = pileups
}

workflow {
    if (!params.samples) {
        error 'Set --samples: a CSV with columns sample, r1, r2'
    }
    if (!params.index && !params.targets) {
        error 'Set --index (a kfp index) or --targets (a FASTA to build one from)'
    }
    def sheet = file(params.samples, checkIfExists: true)
    // reads paths may be relative to the sample sheet
    def reads = channel.fromPath(sheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample, file(sheet.parent.resolve(row.r1), checkIfExists: true),
                file(sheet.parent.resolve(row.r2), checkIfExists: true))
        }
    def index = params.index ? channel.value(file(params.index, checkIfExists: true))
        : INDEX(file(params.targets, checkIfExists: true), params.index_args).index
    def store = params.store ? channel.fromPath(params.store, checkIfExists: true).collect()
        : channel.value([])
    PILEUP(reads, index, store, file("${projectDir}/assets/NO_WITHHOLD"))
}
