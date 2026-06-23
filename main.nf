#!/usr/bin/env nextflow

process BUILD_SAMPLESHEET { 

    container params.metamist_container

    output:
        path "*.tsv", emit: samplesheet

    script:
    def cohort_list = params.cohorts.replace(',', ' ')
    """
    fetch_cohort_samplesheet.py --project ${params.metamist_project} --cohorts ${cohort_list}
    """
    
}

process FILTER_AND_CONVERT {

    container params.bcftools_container

    input:
        tuple val(meta), path(gvcf)
    
    output:
        tuple val(meta), path("${meta.id}.converted.vcf.gz"), emit: vcf
    
    
    script:
    """
    bcftools norm -m -any $gvcf \
        | bcftools view -e 'ALT="<NON_REF>" || FMT/DP="." || FMT/DP < ${params.min_depth} || FMT/GQ < ${params.min_gq}' \
            --trim-alt-alleles -Oz -o ${meta.id}.converted.vcf.gz
    """
}

process RBCEQ2 {

    container params.rbceq2_container

    publishDir { "${params.outdir}/${params.output_version}/${meta.id}" }, mode: 'copy' // after success, copy declared outputs into results/

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("${meta.id}_*.tsv"), emit: results // globs all three .tsv files 
        path("${meta.id}_PDFs"), optional: true, emit: pdfs
    
    script:
    def output_pdfs = params.output_pdfs ? '--PDFs' : ''
    """
    rbceq2 \
        --vcf $vcf \
        --out ${meta.id} \
        --reference_genome ${params.reference_genome} \
        ${output_pdfs}
    """
}

process GATHER {

    publishDir {"${params.outdir}/${params.output_version}/${cohort}/combined"}, mode: 'copy'

    input:
        tuple val(cohort), path(tsvs)

    output:
        tuple val(cohort), path("combined.${cohort}.*.tsv"), emit: combined

    script:
    """
    awk 'NR==1 || FNR>1' *_geno.tsv > combined.${cohort}.geno.tsv
    awk 'NR==1 || FNR>1' *_pheno_numeric.tsv > combined.${cohort}.pheno_numeric.tsv 
    awk 'NR==1 || FNR>1' *_pheno_alphanumeric.tsv > combined.${cohort}.pheno_alphanumeric.tsv 
    """
}

process REGISTER_METAMIST {

    tag "${cohort}"
    container params.metamist_container

    input:
        tuple val(cohort), path(tsvs)

    script:
    def out_prefix = "${params.outdir}/${params.output_version}/${cohort}/combined"
        """
        update_metamist.py \
            --project ${params.metamist_project} \
            --cohorts ${cohort} \
            --type custom \
            --output ${out_prefix}/combined.${cohort}.geno.tsv \
            --secondary pheno_numeric=${out_prefix}/combined.${cohort}.pheno_numeric.tsv pheno_alphanumeric=${out_prefix}/combined.${cohort}.pheno_alphanumeric.tsv

        """
}

workflow {

    BUILD_SAMPLESHEET()
    gvcf_ch = BUILD_SAMPLESHEET.out.samplesheet
        .splitCsv(header: true, sep: '\t')
        .map { row -> [[id: row.sg_id, cohort: row.cohort, project: row.project], file(row.gvcf)] }
    
    FILTER_AND_CONVERT(gvcf_ch)

    RBCEQ2(FILTER_AND_CONVERT.out.vcf)

    cohort_gathered_ch = RBCEQ2.out.results
        .map { meta, files -> [meta.cohort, files] }   // key each item by cohort
        .groupTuple()                                   // fan-in per distinct cohort
        .map { cohort, fileLists -> [cohort, fileLists.flatten()] }  // 3 TSVs × N samples → flat list

    GATHER(cohort_gathered_ch)

    REGISTER_METAMIST(GATHER.out.combined)

}