#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_dir = "trim_fastp_exo"
params.outdir    = "multiqc_report"

workflow {

    // Channel with all fastp HTML files
    html_files = Channel.fromPath("${params.input_dir}/*_fastp.html")

    run_multiqc(html_files.collect())
}

process run_multiqc {

    input:
    path qc_files

    output:
    path "${params.outdir}"

    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    mkdir -p ${params.outdir}
    multiqc ${qc_files.join(' ')} -o ${params.outdir}
    """
}
