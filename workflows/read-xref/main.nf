// Read cross-reference pipeline (plan: Genome profile and functional update, Read
// cross-reference pipeline; phase 12): each read classified twice, independently, then
// joined on read id into a function (MGnify90 unit) x genome table. Per sample:
// read -> function: workflows/pileup (kfp query, recruited reference, DIAMOND or bwa-mem2,
// read posteriors); read -> genome: sylph picks the genomes present, bowtie2 maps every read
// competitively to them; join: python -m kfp_genomes.xref. See README.md.

include { INDEX; PILEUP } from '../pileup/main.nf'

process SKETCH {
    label 'process_medium'
    container 'quay.io/biocontainers/sylph:0.9.0--ha6fb395_0'

    input:
    path fastas

    output:
    path 'genomes.syldb', emit: db

    script:
    "sylph sketch -g ${fastas} -o genomes -t ${task.cpus} ${params.sketch_args}"

    stub:
    "touch genomes.syldb"
}

process SYLPH {
    tag "${sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/sylph:0.9.0--ha6fb395_0'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/${it}" }

    input:
    tuple val(sample), path(r1), path(r2)
    path db

    output:
    tuple val(sample), path('sylph.tsv'), emit: profile

    script:
    "sylph profile ${db} -1 ${r1} -2 ${r2} -t ${task.cpus} ${params.sylph_args} > sylph.tsv"

    stub:
    "touch sylph.tsv"
}

process BOWTIE2 {
    tag "${sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/bowtie2:2.5.5--h63e9258_1'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/${it}" },
        pattern: 'contigs.tsv'

    input:
    tuple val(sample), path(r1), path(r2), path(sylph)
    path fastas
    path names  // file name -> genome

    output:
    tuple val(sample), path('hits.tsv'), path('contigs.tsv'), emit: hits

    script:
    // The genomes sylph reports (or every genome with --sylph false), concatenated with
    // contigs renamed genome|n; then hits as TSV (kfp_genomes.xref.HITS_COLUMNS).
    def pick = params.sylph
        ? "tail -n +2 ${sylph} | cut -f2 | xargs -r -n1 basename"
        : "cut -f1 ${names}"
    """
    ${pick} | sort -u > picked.txt
    awk -F'\\t' 'NR==FNR { want[\$1]; next } \$1 in want' picked.txt ${names} |
    while IFS=\$'\\t' read -r f g; do
        case \$f in *.gz) gzip -dc \$f ;; *) cat \$f ;; esac |
            awk -v g="\$g" '/^>/ { print ">" g "|" n++; next } { print }'
    done > ref_genomes.fa
    awk 'BEGIN { OFS = "\\t"; print "contig", "genome", "length" }
        /^>/ { if (c) print c, g, l; c = substr(\$0, 2); g = c; sub(/\\|[^|]*\$/, "", g); l = 0; next }
        { l += length(\$0) } END { if (c) print c, g, l }' ref_genomes.fa > contigs.tsv
    touch hits.tsv
    if [ -s ref_genomes.fa ]; then
        bowtie2-build --threads ${task.cpus} ref_genomes.fa ref_genomes > /dev/null
        bowtie2 -x ref_genomes -1 ${r1} -2 ${r2} -k ${params.bowtie2_k} --no-unal \\
            -p ${task.cpus} ${params.bowtie2_args} |
        awk 'BEGIN { OFS = "\\t" } /^@/ { next } {
            a = 0; c = \$6
            while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
                n = substr(c, 1, RLENGTH - 1) + 0; op = substr(c, RLENGTH, 1)
                if (op ~ /[MI=X]/) a += n
                c = substr(c, RLENGTH + 1)
            }
            s = 0; m = 0
            for (i = 12; i <= NF; i++) {
                if (\$i ~ /^AS:i:/) s = substr(\$i, 6)
                if (\$i ~ /^NM:i:/) m = substr(\$i, 6)
            }
            print \$1, \$2, \$3, \$5, s, m, a
        }' > hits.tsv
    fi
    """

    stub:
    "touch hits.tsv contigs.tsv"
}

process XREF {
    tag "${sample}"
    label 'process_low'
    publishDir params.outdir, mode: 'copy', saveAs: { "${sample}/${it}" }

    input:
    tuple val(sample), path(pileup), path(ref), path(hits), path(contigs)
    path genomes

    output:
    tuple val(sample), path('xref'), emit: xref

    script:
    """
    ${params.venv}/bin/python -m kfp_genomes.xref ${pileup} ${ref} ${hits} ${contigs} xref \\
        --genomes ${genomes} --min-identity ${params.min_identity} --min-mapq ${params.min_mapq}
    """

    stub:
    "mkdir xref"
}

workflow {
    if (!params.samples || !params.genomes) {
        error 'Set --samples (CSV: sample, r1, r2) and --genomes (TSV: genome, path, optional taxonomy)'
    }
    if (!params.index && !params.targets) {
        error 'Set --index (a kfp index) or --targets (a FASTA to build one from)'
    }
    def sheet = file(params.samples, checkIfExists: true)
    def reads = channel.fromPath(sheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample, file(sheet.parent.resolve(row.r1), checkIfExists: true),
                file(sheet.parent.resolve(row.r2), checkIfExists: true))
        }
    def genome_sheet = file(params.genomes, checkIfExists: true)
    def genomes = channel.fromPath(genome_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.genome, file(genome_sheet.parent.resolve(row.path), checkIfExists: true)) }
    // staged by file name, so file names must be unique
    def fastas = genomes.map { genome, fasta -> fasta }.collect()
    def names = genomes.map { genome, fasta -> "${fasta.name}\t${genome}" }
        .collectFile(name: 'genome_names.tsv', newLine: true)
    def index = params.index ? channel.value(file(params.index, checkIfExists: true))
        : INDEX(file(params.targets, checkIfExists: true)).index
    def store = params.store ? channel.fromPath(params.store, checkIfExists: true).collect()
        : channel.value([])

    def pileup = PILEUP(reads, index, store)
    def selected = params.sylph
        ? SYLPH(reads, params.sylph_db ? file(params.sylph_db, checkIfExists: true)
            : SKETCH(fastas).db).profile
        : reads.map { sample, r1, r2 -> tuple(sample, genome_sheet) }  // a placeholder, unused
    def hits = BOWTIE2(reads.join(selected), fastas, names).hits
    XREF(pileup.pileup.join(pileup.ref).join(hits), genome_sheet)
}
