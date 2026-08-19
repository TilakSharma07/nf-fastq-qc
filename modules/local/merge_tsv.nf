process MERGE_TSV {
    tag "${name}"
    label 'process_single'
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    tuple val(name), path(tsvs)

    output:
    path "${name}.tsv", emit: tsv

    script:
    """
    # Header from the first file, then every data row, sorted by sample so the
    # merged table is deterministic regardless of task completion order.
    for f in \$(ls ${tsvs} | sort); do head -n 1 "\$f"; break; done > ${name}.tsv
    for f in \$(ls ${tsvs} | sort); do tail -n +2 "\$f"; done | sort -k1,1 >> ${name}.tsv
    """
}
