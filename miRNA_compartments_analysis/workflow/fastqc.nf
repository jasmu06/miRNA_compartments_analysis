#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = 'trimmed_reads/tissue_trimmed/*.fq.gz'

params.outdir = 'fastqc_reports'

process run_fastqc {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path fastq_file

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc ${fastq_file}
    """
}

workflow {
    Channel
        .fromPath(params.input)
        .set { fastq_files }

    run_fastqc(fastq_files)
}
