#!/usr/bin/env nextflow

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'
include { FILTER_AND_CONVERT } from './modules/local/filter_and_convert'
include { RBCEQ2 }            from './modules/local/rbceq2'
include { GATHER }            from './modules/local/gather'
include { REGISTER_METAMIST } from './modules/local/register_metamist'

workflow {

    BUILD_SAMPLESHEET()
    gvcf_ch = BUILD_SAMPLESHEET.out.samplesheet
        .splitCsv(header: true, sep: '\t')
        .map { row -> [[id: row.sg_id, cohort: row.cohort, project: row.project], file(row.gvcf)] }
    
    FILTER_AND_CONVERT(gvcf_ch, params.min_depth, params.min_gq)

    RBCEQ2(FILTER_AND_CONVERT.out.vcf)

    cohort_gathered_ch = RBCEQ2.out.results
        .map { meta, files -> [meta.cohort, files] }   // key each item by cohort
        .groupTuple()                                   // fan-in per distinct cohort
        .map { cohort, fileLists -> [cohort, fileLists.flatten()] }  // 3 TSVs × N samples → flat list

    GATHER(cohort_gathered_ch)

    REGISTER_METAMIST(GATHER.out.combined)

}