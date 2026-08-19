process FASTQ_TRIM {
    tag "${meta.id}"
    label 'process_single'
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*.trimmed.fastq.gz"
    publishDir "${params.outdir}/trim_logs", mode: 'copy', pattern: "*.trim.log.tsv"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.trim.log.tsv")    , emit: log
    path "versions.yml"                                 , emit: versions

    script:
    """
    fq_trim.py \\
        --fastq ${reads} \\
        --sample ${meta.id} \\
        --out-fastq ${meta.id}.trimmed.fastq.gz \\
        --out-log ${meta.id}.trim.log.tsv \\
        --quality-cutoff ${params.quality_cutoff} \\
        --min-length ${params.min_length} \\
        --adapter ${params.adapter}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
