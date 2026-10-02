#!/usr/bin/env nextflow

// This process runs Bambu on the sorted and indexed BAM files produced by the previous steps of the workflow. 
// The output will be a directory named "bambu_output" containing the Bambu results for each barcode and condition.


 process runBambu {
    tag "Bambu_quantification"
    conda "${params.env_dir}/bambu.yml"
    cpus params.bambu_threads

    publishDir "${params.outdir}/bambu", mode: 'copy'

    input:
    path bams
    path reference_genome
    path annotation_gtf

    output:
    path "bambu_output/*"

    script:

    """
    mkdir -p bambu_output

    Rscript - <<'RSCRIPT'

    library(bambu)

    message("========================================")
    message("Bambu quantification")
    message("========================================")

    annotations <- prepareAnnotations(
        "${annotation_gtf}"
    )

    bam_files <- c(
        ${bams.collect { "\"${it}\"" }.join(",\n        ")}
    )

    message("BAM files:")
    message(paste(bam_files, collapse = "\\n"))

    message("GTF: ${annotation_gtf}")
    message("Genome: ${reference_genome}")

    message("Discovery: FALSE")
    message("Quantification: TRUE")
    message("Low memory: ${params.bambu_low_memory}")
    message("Threads: ${task.cpus}")

    message("========================================")

    se <- bambu(
        reads = bam_files,
        annotations = annotations,
        genome = "${reference_genome}",
        discovery = FALSE,
        quant = TRUE,
        lowMemory = ${params.bambu_low_memory ? 'TRUE' : 'FALSE'},
        ncore = ${task.cpus}
    )

    writeBambuOutput(
        se,
        path = "bambu_output"
    )

    message("========================================")
    message("Bambu completed successfully")
    message("========================================")

    RSCRIPT
    """
}