#!/usr/bin/env nextflow

// Converting SAM files to BAM files and sorting.


process runSamtoBam {
    tag "${barcode}_${condition}"
    conda "${params.env_dir}/samtools.yml"

    input:
    tuple val(barcode), val(condition), path(sam)
    output:
    tuple val(barcode), val(condition),
          path("barcode_${barcode}_${condition}_aln.bam")
    script:
    """
    echo "Converting ${barcode} (${condition}) SAM to BAM"
    samtools view \
        -b \
        ${sam} \
        > barcode_${barcode}_${condition}_aln.bam
    """
}


process runSortBam {
    conda "${params.env_dir}/samtools.yml"

    input:
    tuple val(barcode), val(condition), path(bam)
    output:
    tuple val(barcode), val(condition), path("barcode_${barcode}_${condition}_aln_sorted.bam")

    script:
    """
    echo "Starting Samtools conversion"
    echo "Converting $barcode SAM to BAM:"
    samtools sort $bam -o barcode_${barcode}_${condition}_aln_sorted.bam
    """
}


process runBamIndex {
    tag "${barcode}_${condition}"
    conda "${params.env_dir}/samtools.yml"
    publishDir "${params.outdir}/bams",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            "${barcode}_${condition}/${filename}"
        }

    input:
    tuple val(barcode), val(condition), path(bam)
    output:
    tuple val(barcode), val(condition),
          path(bam),
          path("${bam}.bai")

    script:
    """
    echo "Indexing ${barcode} (${condition}) BAM"
    samtools index ${bam}
    """
}
