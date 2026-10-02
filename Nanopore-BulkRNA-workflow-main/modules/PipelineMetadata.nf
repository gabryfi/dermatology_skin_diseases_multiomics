#!/usr/bin/env nextflow

process savePipelineParameters {

    publishDir "${params.outdir}/pipeline_info",
        mode: 'copy',
        overwrite: true

    output:
    path 'parameters.json'

    script:
    def parameters = groovy.json.JsonOutput.prettyPrint(
        groovy.json.JsonOutput.toJson(params)
    )

    """
    cat > parameters.json <<'EOF'
${parameters}
EOF
    """
}

// WIP quando avroò tutti gli ambienti va inserito, ora printa solo nextflow version
process saveSoftwareVersions {

    publishDir "${params.outdir}/pipeline_info",
        mode: 'copy',
        overwrite: true

    output:
    path 'software_versions.txt'

    script:
    """
    {
        echo "========================================"
        echo "Pipeline software versions"
        echo "========================================"

        echo ""
        echo "Nextflow:"
        nextflow -version 2>&1 | head -n 5

    } > software_versions.txt
    """
}