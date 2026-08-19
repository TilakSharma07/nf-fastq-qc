process VALIDATE_QC {
    label 'process_single'
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path raw_qc
    path trimmed_qc
    path trim_log

    output:
    path "validation.txt", emit: report

    script:
    """
    validate_qc.py \\
        --raw-qc ${raw_qc} \\
        --trimmed-qc ${trimmed_qc} \\
        --trim-log ${trim_log} \\
        --min-q30-after ${params.min_q30_after} \\
        --out validation.txt
    """
}
