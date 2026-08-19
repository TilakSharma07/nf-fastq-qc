process FASTQ_QC {
    tag "${meta.id}:${meta.stage}"
    label 'process_single'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.${meta.stage}.qc.tsv") , emit: tsv
    tuple val(meta), path("${meta.id}.${meta.stage}.qc.json"), emit: json
    path "versions.yml"                                      , emit: versions

    script:
    """
    fq_qc.py \\
        --fastq ${reads} \\
        --sample ${meta.id} \\
        --out-tsv ${meta.id}.${meta.stage}.qc.tsv \\
        --out-json ${meta.id}.${meta.stage}.qc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
