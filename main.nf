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

process CONVERT {

    container params.bcftools_container

    input:
        tuple val(meta), path(gvcf)
    
    output:
        tuple val(meta), path("${meta.id}.converted.vcf.gz"), emit: vcf
    
    
    script:
    """
    bcftools norm -m -any $gvcf | bcftools view -e 'ALT="<NON_REF>"' --trim-alt-alleles -Oz -o ${meta.id}.converted.vcf.gz
    """
}

process RBCEQ2 {

    container params.rbceq2_container

    publishDir { "${params.outdir}/${meta.id}" }, mode: 'copy' // after success, copy declared outputs into results/

    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("${meta.id}_*.tsv"), emit: results // globs all three .tsv files 
    
    script:
    """
    rbceq2 --vcf $vcf --out ${meta.id} --reference_genome GRCh38
    """
}

workflow {

    BUILD_SAMPLESHEET()
    gvcf_ch = BUILD_SAMPLESHEET.out.samplesheet
        .splitCsv(header: true, sep: '\t')
        .map { row -> [[id: row.sg_id, cohort: row.cohort, project: row.project], file(row.gvcf)] }
    
    CONVERT(gvcf_ch)

    RBCEQ2(CONVERT.out.vcf)
}