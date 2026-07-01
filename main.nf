#!/usr/bin/env nextflow

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'
include { FILTER_AND_CONVERT } from './modules/local/filter_and_convert'
include { RBCEQ2 }            from './modules/local/rbceq2'
include { GATHER }            from './modules/local/gather'
include { REGISTER_METAMIST } from './modules/local/register_metamist'

workflow {

    BUILD_SAMPLESHEET(params.metamist_project, params.cohorts)

    gvcf_ch = BUILD_SAMPLESHEET.out.samplesheet
        .splitCsv(header: true, sep: '\t')
        .map { row -> [[id: row.sg_id, cohort: row.cohort, project: row.project], file(row.gvcf), file(row.gvcf + '.tbi')] }

    FILTER_AND_CONVERT(gvcf_ch, params.region_bed, params.min_depth, params.min_gq, params.n_cpus)

    RBCEQ2(FILTER_AND_CONVERT.out.vcf, params.reference_genome, params.output_pdfs)

    cohort_gathered_ch = RBCEQ2.out.results
        .map { meta, files -> [meta.cohort, files] }   // key each item by cohort
        .groupTuple()                                   // fan-in per distinct cohort
        .map { cohort, fileLists -> [cohort, fileLists.flatten()] }  // 3 TSVs × N samples → flat list

    GATHER(cohort_gathered_ch)

    REGISTER_METAMIST(GATHER.out.combined, params.metamist_project, params.outdir, params.output_version)

}