process MAKE_REPORT {
    label 'process_single'
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path raw_qc
    path trimmed_qc
    path trim_log
    path validation
    path figure

    output:
    path "qc_report.md", emit: report

    script:
    """
    make_report.py \\
        --raw-qc ${raw_qc} \\
        --trimmed-qc ${trimmed_qc} \\
        --trim-log ${trim_log} \\
        --validation ${validation} \\
        --figure ${figure} \\
        --out qc_report.md
    """
}
