#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQ_QC as QC_RAW     } from './modules/local/fastq_qc'
include { FASTQ_QC as QC_TRIMMED } from './modules/local/fastq_qc'
include { FASTQ_TRIM             } from './modules/local/fastq_trim'
include { MERGE_TSV as MERGE_RAW } from './modules/local/merge_tsv'
include { MERGE_TSV as MERGE_TRIM} from './modules/local/merge_tsv'
include { MERGE_TSV as MERGE_LOG } from './modules/local/merge_tsv'
include { VALIDATE_QC            } from './modules/local/validate_qc'
include { PLOT_QC                } from './modules/local/plot_qc'
include { MAKE_REPORT            } from './modules/local/make_report'

def helpMessage() {
    log.info """
    nf-fastq-qc — reproducible FASTQ QC and trimming with a hard validation gate

    Usage:
      nextflow run . --input samplesheet.csv --outdir results -profile docker

    Required:
      --input               CSV samplesheet with columns: sample,fastq
      --outdir              Output directory

    Trimming options:
      --quality-cutoff      Phred cutoff for 3' trimming     [${params.quality_cutoff}]
      --min_length          Discard reads shorter than this   [${params.min_length}]
      --adapter             Adapter sequence to remove        [${params.adapter}]

    Validation:
      --min_q30_after       Required %Q30 floor after trimming [${params.min_q30_after}]
    """.stripIndent()
}

workflow {

    if (params.help) {
        helpMessage()
        System.exit(0)
    }

    if (!params.input) {
        error "No samplesheet provided. Use --input samplesheet.csv (see --help)."
    }

    // Parse the samplesheet into (meta, reads) tuples, failing loudly on a
    // missing FASTQ rather than silently dropping the sample.
    ch_reads = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample?.trim()) {
                error "Samplesheet row is missing a 'sample' value: ${row}"
            }
            if (!row.fastq?.trim()) {
                error "Sample '${row.sample}' is missing a 'fastq' path."
            }
            def fq = file(row.fastq, checkIfExists: true)
            tuple([id: row.sample.trim()], fq)
        }

    QC_RAW(ch_reads.map { meta, reads -> tuple(meta + [stage: 'raw'], reads) })

    FASTQ_TRIM(ch_reads)

    QC_TRIMMED(
        FASTQ_TRIM.out.reads.map { meta, reads -> tuple(meta + [stage: 'trimmed'], reads) }
    )

    // Collect per-sample tables into single merged summaries.
    MERGE_RAW (QC_RAW.out.tsv.map     { _meta, tsv -> tsv }.collect().map { tuple('raw_qc', it) })
    MERGE_TRIM(QC_TRIMMED.out.tsv.map { _meta, tsv -> tsv }.collect().map { tuple('trimmed_qc', it) })
    MERGE_LOG (FASTQ_TRIM.out.log.map { _meta, lg  -> lg  }.collect().map { tuple('trim_log', it) })

    // Hard gate: a failing invariant stops the run before any report exists.
    VALIDATE_QC(MERGE_RAW.out.tsv, MERGE_TRIM.out.tsv, MERGE_LOG.out.tsv)

    ch_json = QC_RAW.out.json.map { _meta, j -> j }
        .mix(QC_TRIMMED.out.json.map { _meta, j -> j })
        .collect()

    PLOT_QC(MERGE_RAW.out.tsv, MERGE_TRIM.out.tsv, ch_json)

    MAKE_REPORT(
        MERGE_RAW.out.tsv,
        MERGE_TRIM.out.tsv,
        MERGE_LOG.out.tsv,
        VALIDATE_QC.out.report,
        PLOT_QC.out.figure
    )
}
