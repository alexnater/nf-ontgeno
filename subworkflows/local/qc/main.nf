/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../../../modules/nf-core/fastqc'
include { NANOPLOT               } from '../../../modules/nf-core/nanoplot'
include { TOULLIGQC              } from '../../../modules/nf-core/toulligqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QC {

    take:
    ch_reads       // channel: [ meta, fastq ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run NanoPlot
    //
    NANOPLOT (
        ch_reads
    )
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())

    //
    // MODULE: Run ToulligQC
    //
    TOULLIGQC (
        ch_reads
    )
    ch_versions = ch_versions.mix(TOULLIGQC.out.versions.first())

    emit:
    fastqc = FASTQC.out.zip       // channel: [ meta, html ]
    nanoplot = NANOPLOT.out.txt   // channel: [ meta, txt ]
    versions = ch_versions        // channel: [ path(versions.yml) ]
}