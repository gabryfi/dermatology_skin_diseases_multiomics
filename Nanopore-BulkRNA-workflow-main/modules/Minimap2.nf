#!/usr/bin/env nextflow

// Aligning Post QC-ready fastq files.


process runMinimap2 {
    tag "${barcode}_${condition}"
    conda "${params.env_dir}/mapping.yml"

    input:
    tuple val(barcode), val(condition), path(fastq)
    output:
    tuple val(barcode), val(condition),
          path("barcode_${barcode}_${condition}_aln.sam")

    script:
    """
    echo "Starting Minimap2 alignment"
    echo "Aligning ${barcode} (${condition})"
    minimap2 \
        -t ${params.default_threads} \
        -a \
        -x splice \
        ${params.reference_genome} \
        ${fastq} \
        > barcode_${barcode}_${condition}_aln.sam
    """
}