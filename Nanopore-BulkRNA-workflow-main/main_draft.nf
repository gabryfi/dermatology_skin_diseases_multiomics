// --------------------------------------------- //
// Nextflow workflow for Bulk RNA-seq analysis using Nanopore data
// --------------------------------------------- //

// Pipeline parameters
// The parameters are defined in the nextflow.config file, which is located in the root directory of the workflow.

// --------------------------------------------- //
// Processes
// --------------------------------------------- //


// samplesheet example 
// sample,barcode,condition
//  1,barcode_01,veh
//  2,barcode_02,veh
//  3,barcode_03,veh
//  ecc.

// --------------------------------------------- //
// Workflow
// --------------------------------------------- //

include { validateReferences } from './modules/ValidateReferences.nf'
include { saveSoftwareVersions } from './modules/PipelineMetadata.nf'
include { savePipelineParameters } from './modules/PipelineMetadata.nf'
include { saveSamplesheet } from './modules/SaveSamplesheet.nf'
include { readSamplesheet } from './modules/ReadSamplesheet.nf'
include { runConcatenateFastq } from './modules/ConcatenateFastq.nf'
include { runNanoPlotQC_pre } from './modules/NanoPlotQC.nf'
include { runNanoPlotQC_post } from './modules/NanoPlotQC.nf'
include { aggregateNanoPlotQC } from './modules/aggregateNanoPlot.nf'
include { runFastCat } from './modules/FastCat_filtering.nf'  
include { runPychopper } from './modules/Pychopper.nf'  
include { runMinimap2 } from './modules/Minimap2.nf'
include { runSamtoBam } from './modules/Samtools.nf'
include { runSortBam } from './modules/Samtools.nf' 
include { runBamIndex } from './modules/Samtools.nf'    
include { runBambu } from './modules/Bambu.nf'

workflow {
    samplesheet_ch = channel.fromPath(params.samplesheet)
    saveSamplesheet(samplesheet_ch)
    readSamplesheet(samplesheet_ch)
    samplesheet_ch
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.barcode,
                row.condition,
                file("${params.fastq_dir}/${row.barcode}/*.fastq.gz", checkIfExists: true)
            )
        }
        .view { barcode, condition, fastq_files -> "INPUT: ${barcode} | condition: ${condition} | FASTQ files: ${fastq_files.size()}"
        }
        .set { fastq_ch }   
    genome_ch = channel.fromPath(params.reference_genome, checkIfExists: true)
    annotation_ch = channel.fromPath(params.annotation_file, checkIfExists: true)
    checksum_ch = channel.fromPath("${baseDir}/references/checksum.sha256", checkIfExists: true)

    validateReferences(genome_ch, annotation_ch, checksum_ch)
    runConcatenateFastq(fastq_ch)
    runNanoPlotQC_pre(runConcatenateFastq.out)
    runFastCat(runConcatenateFastq.out)
    runNanoPlotQC_post(runFastCat.out.filtered)
    aggregateNanoPlotQC(runNanoPlotQC_pre.out.collect(), runNanoPlotQC_post.out.collect())
    runPychopper(runFastCat.out.filtered)
    runMinimap2(runPychopper.out.full_length)
    runSamtoBam(runMinimap2.out)
    runSortBam(runSamtoBam.out)
    runBamIndex(runSortBam.out)
	bam_ch = runBamIndex.out.map { barcode, condition, bam, bai -> bam }.collect()
    runBambu(bam_ch, file(params.reference_genome), file(params.annotation_file))
    savePipelineParameters()
    saveSoftwareVersions()
}

// CLI final command to run the workflow
// nextflow run main.nf --samplesheet samplesheet.csv --fastq_dir /data/fastq_pass --reference_genome 









