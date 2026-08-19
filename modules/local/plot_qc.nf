process PLOT_QC {
    label 'process_single'
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path raw_qc
    path trimmed_qc
    path json_files

    output:
    path "qc_summary.png", emit: figure

    script:
    """
    mkdir -p jsons
    cp ${json_files} jsons/
    plot_qc.py \\
        --raw-qc ${raw_qc} \\
        --trimmed-qc ${trimmed_qc} \\
        --json-dir jsons \\
        --q30-floor ${params.min_q30_after} \\
        --out qc_summary.png
    """
}
