/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                  } from '../../../modules/nf-core/fastp'
include { FASTPLONG              } from '../../../modules/nf-core/fastplong'
include { CHOPPER                } from '../../../modules/nf-core/chopper'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPROCESSING {

    take:
    ch_reads       // channel: [ meta, fastq ]
    trimmer        // string

    main:

    ch_versions = Channel.empty()
    ch_json = Channel.empty()

    if (trimmer == 'fastp') {
        //
        // MODULE: Run fastp
        //
        FASTP (
            ch_reads,
            [],
            false,
            false,
            false
        ).reads
         .set { trimmed }
        ch_json = ch_json.mix(FASTP.out.json)
        ch_versions = ch_versions.mix(FASTP.out.versions.first())

    } else if (trimmer == 'fastplong') {
        //
        // MODULE: Run fastplong
        //
        FASTPLONG (
            ch_reads,
            [],
            false,
            false
        ).reads
         .set { trimmed }
        ch_json = ch_json.mix(FASTPLONG.out.json)
        ch_versions = ch_versions.mix(FASTPLONG.out.versions.first())

    } else {
        //
        // MODULE: Run chopper
        //
        CHOPPER (
            ch_reads,
            []
        ).fastq
         .set { trimmed }
        ch_versions = ch_versions.mix(CHOPPER.out.versions.first())
    }

    emit:
    trimmed                       // channel: [ meta, fastq ]
    json = ch_json                // channel: [ meta, json ]
    versions = ch_versions        // channel: [ path(versions.yml) ]
}