#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "fastq_exo/*.fastq.gz"


params.outdir = "trim_fastp_exo"

workflow {
    Channel
        .fromPath(params.input)
        .map { file -> tuple(file.baseName, file) }
        .set { raw_reads }

    trim_fastp(raw_reads)
}

process trim_fastp {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    path("${sample_id}_trimmed.fastq"), emit: trimmed
    path("${sample_id}_fastp.html"), emit: html
    path("${sample_id}_fastp.json"), emit: json

    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    fastp -i ${fastq_file} \\
          -o ${sample_id}_trimmed.fastq \\
          --length_required 18 \\
          --length_limit 30 \\
          --qualified_quality_phred 20 \\
          --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \\
          --html ${sample_id}_fastp.html \\
          --json ${sample_id}_fastp.json
    """
}
