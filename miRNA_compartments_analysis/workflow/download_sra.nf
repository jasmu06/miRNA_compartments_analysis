#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.srr_list = '/Users/jasmu/mirna_analysis/exomissed.txt'

workflow {
    Channel
        .fromPath(params.srr_list)
        .splitText()
        .map { it.trim() }
        .filter { it }
        .set { srr_ids }

    print_srr(srr_ids)
    download_fastq(srr_ids)
}

process print_srr {
    input:
    val srr

    output:
    stdout()

    script:
    """
    echo "SRR to download: $srr"
    """
}

process download_fastq {
    input:
    val srr

    output:
    file("${srr}*")

    script:
    """
    echo "Downloading $srr"
    fasterq-dump $srr --outdir . --temp ./tmp
    echo "Files created in:"
    ls -lh
    """
}
