// aai-rescore: a finished fmh-benchmark run's MGnify unit profiles rescored against
// nearest-member identity, with aai re-estimated under survival models (plan, phase 7,
// step 34), without rebuilding its indexes or rerunning its queries. Steps are bench.py's.

include { MGNIFY_DB; MGNIFY_ANNOTATE } from '../fmh-benchmark/mgnify_diamond.nf'

process PROTEINS {
    label 'process_single'

    input:
    path genomes

    output:
    path 'proteins.*.faa'

    script:
    "${params.bench} pfam-proteins --genomes-dir ${genomes} --chunks ${params.chunks}"

    stub:
    "touch proteins.0.faa proteins.1.faa"
}

process MEMBERS {
    label 'process_medium'

    input:
    path gene_units
    path members, stageAs: 'members/*'

    output:
    tuple val('members'), path('members.faa'), emit: faa
    path 'member_clusters.parquet', emit: clusters

    script:
    "${params.bench} mgnify-members --members ${members} --gene-units ${gene_units}"

    stub:
    "touch members.faa member_clusters.parquet"
}

process NEAREST {
    label 'process_medium'
    publishDir params.outdir, mode: 'copy'

    input:
    path hits, stageAs: 'hits/*'
    path gene_units, stageAs: 'reps.parquet'
    path clusters

    output:
    path 'gene_units.parquet'

    script:
    """
    ${params.bench} mgnify-nearest --hits hits/* --gene-units reps.parquet \\
        --member-clusters ${clusters} --out gene_units.parquet
    """

    stub:
    "touch gene_units.parquet"
}

process SCORE {
    tag "seed ${sid} ${index}${arm ? '~' + arm : ''}${model_name ? '+' + model_name : ''}"
    label 'process_single'

    input:
    tuple val(sid), val(index), val(arm), path(units), path(genes), val(model_name), val(model)
    path gene_units

    output:
    path 'aai_score.tsv'

    script:
    def scored = model_name ? "${arm}+${model_name}" : arm
    """
    ${params.bench} aai-score --profile ${units} --gene-units ${gene_units} --genes ${genes} \\
        --sample seed${sid} --index ${index} --arm '${scored}' ${model ? "--model ${model}" : ''}
    """

    stub:
    "touch aai_score.tsv"
}

process SUMMARY {
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

workflow {
    def run = file(params.run, checkIfExists: true)
    def reps = file(params.gene_units ?: "${run}/mgnify/${params.annotation}_gene_units.parquet", checkIfExists: true)
    MEMBERS(reps, files(params.members, checkIfExists: true))
    MGNIFY_ANNOTATE(
        MGNIFY_DB(MEMBERS.out.faa).combine(PROTEINS(file(params.genomes, checkIfExists: true)).flatten()),
        params.nearest_diamond_args,
    )
    NEAREST(MGNIFY_ANNOTATE.out.hits.map { _name, hits -> hits }.collect(), reps, MEMBERS.out.clusters)
    // units/seed<sid>_<index>[~<arm>].tsv, as the fmh-benchmark run published them
    def ch_units = channel.fromPath("${run}/units/seed*.tsv").map { f ->
        def m = (f.baseName =~ /^seed([^_]+)_([^~]+)(?:~(.*))?$/)
        if (!m.matches()) {
            error "Unexpected unit profile name: ${f.name}"
        }
        [m.group(1), m.group(2), m.group(3) ?: '', f, file("${run}/truth/seed${m.group(1)}_genes.csv", checkIfExists: true)]
    }
    def ch_models = channel.of(['', '']).mix(
        channel.fromList(params.models).map { m -> [m.name, file(m.path, checkIfExists: true).toString()] }
    )
    SUMMARY(SCORE(ch_units.combine(ch_models), NEAREST.out).collect())
}
