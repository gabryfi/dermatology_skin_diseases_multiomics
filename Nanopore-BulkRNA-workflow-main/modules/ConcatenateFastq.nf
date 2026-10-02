#!/usr/bin/env nextflow

// This process, after reading the samplesheet, will concatenate all fastq files for each barcode and condition into a single fastq file. 
// The output will be named as "${barcode}_${condition}_combined.fastq".
// Assuming is requiried to unzip files before concatenation, the process will use zcat to concatenate the fastq files.

process runConcatenateFastq {
    tag "${barcode}_${condition}"

    input:
    tuple val(barcode), val(condition), path(fastq_files)
    output:
    tuple val(barcode), val(condition),
          path("${barcode}_${condition}_combined.fastq")

    script:
    """
    zcat ${fastq_files} > ${barcode}_${condition}_combined.fastq
    """
}
