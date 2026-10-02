#!/usr/bin/env nextflow

process aggregateNanoPlotQC {

    tag "All_samples"
    conda "${params.env_dir}/nanoplot_summary.yml"
    publishDir "${params.outdir}/nanoplot_summary",
        mode: 'copy',
        overwrite: true

    input:
    path pre_dirs
    path post_dirs

    output:
    path "nanoplot_summary"

    script:
    """
    mkdir -p nanoplot_summary

    python3 ${params.env_dir}/../bin/aggregate_nanoplot.py \
        --pre ${pre_dirs} \
        --post ${post_dirs} \
        --outdir nanoplot_summary
    """
}
