#!/usr/bin/env nextflow

process readSamplesheet {
    input:
    path samplesheet
    output:
    path 'samplesheet_check.txt'

    script:
    """
    echo "Samplesheet:"
    cat ${samplesheet} > samplesheet_check.txt
    """
}