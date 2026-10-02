#!/usr/bin/env nextflow

// This process, after loading concatenated fastq files, will perform Filtering QC with Fastcat.


process runFastCat {

    tag "${barcode}_${condition}"
    conda "${params.env_dir}/fastcat.yml"

    input:
    tuple val(barcode), val(condition), path(fastq)
    output:
    tuple val(barcode), val(condition),
          path("${barcode}_${condition}_filtered.fastq"),
          emit: filtered

    script:
    """
    fastcat fastq \
        --min-qscore 10 \
        --min-length 200 \
        -o fastcat_qc \
        ${fastq} \
        > ${barcode}_${condition}_filtered.fastq
    """
}