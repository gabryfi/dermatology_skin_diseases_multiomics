process saveSamplesheet {

    publishDir "${params.outdir}/pipeline_info",
        mode: 'copy',
        overwrite: true

    input:
    path samplesheet

    output:
    path 'samplesheet_used.csv'

    script:
    """
    cp ${samplesheet} samplesheet_used.csv
    """
}