#!/usr/bin/env nextflow

// This process, after loading concatenated fastq files, will output the QC metrics using NanoPlot for each barcode and plot the results in a separate directory. 
// The output will be named as "${barcode}_${condition}_NanoPlotQC".
// I'm working on providing a single summary report for all barcodes and conditions, but for now, the output will be separate for each barcode and condition.

process runNanoPlotQC_pre {
    tag "${barcode}_${condition}"
    conda "${params.env_dir}/nanoplot.yml"
    
    publishDir "${params.outdir}/nanoplot_pre",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            "${barcode}_${condition}/${filename}"
        }

    input:
    tuple val(barcode), val(condition), path(fastq)
    output:
    path "Nanoplot_${barcode}_${condition}"

    script:
    """
    NanoPlot \
        -t ${params.default_threads} \
        --fastq ${fastq} \
        --outdir Nanoplot_${barcode}_${condition}
    """
}


process runNanoPlotQC_post {
    tag "${barcode}_${condition}"
    conda "${params.env_dir}/nanoplot.yml"
    publishDir "${params.outdir}/nanoplot_post",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            "${barcode}_${condition}/${filename}"
        }

    input:
    tuple val(barcode), val(condition), path(fastq)
    output:
    path "Nanoplot_${barcode}_${condition}_filtered"

    script:
    """
    NanoPlot \
        -t ${params.default_threads} \
        --fastq ${fastq} \
        --outdir Nanoplot_${barcode}_${condition}_filtered
    """
}