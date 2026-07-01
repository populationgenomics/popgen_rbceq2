#!/usr/bin/env nextflow

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'
include { FILTER_AND_CONVERT } from './modules/local/filter_and_convert'
include { RBCEQ2 }            from './modules/local/rbceq2'
include { GATHER }            from './modules/local/gather'
include { REGISTER_METAMIST } from './modules/local/register_metamist'

workflow {

    if (params.use_metamist) {
        def cohort_list = params.cohorts.tokenize(',')
        BUILD_SAMPLESHEET(params.metamist_project, cohort_list)
        gvcf_ch = BUILD_SAMPLESHEET.out.samplesheet
            .splitCsv(header: true, sep: '\t')
            .map { row -> [[id: row.sg_id, cohort: row.cohort, project: row.project], file(row.gvcf), file(row.gvcf + '.tbi')] }
    } else {
        gvcf_ch = Channel.fromPath(params.gvcfs)
            .map { f -> [[id: f.simpleName, cohort: params.local_cohort, project: 'local'], f, file(f + '.tbi')] }
    }

    FILTER_AND_CONVERT(gvcf_ch, params.region_bed, params.min_depth, params.min_gq)
    RBCEQ2(FILTER_AND_CONVERT.out.vcf, params.reference_genome)

    cohort_gathered_ch = RBCEQ2.out.results
        .map { meta, files -> [meta.cohort, files] }
        .groupTuple()
        .map { cohort, fileLists -> [cohort, fileLists.flatten()] }

    GATHER(cohort_gathered_ch)

    if (params.use_metamist) {
        REGISTER_METAMIST(GATHER.out.combined, params.metamist_project, params.outdir, params.output_version)
    }
}