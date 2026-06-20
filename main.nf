#!/usr/bin/env nextflow

params.gvcfs = 'test_data/gvcfs/*.g.vcf.gz'

process CONVERT {

    container 'staphb/bcftools:1.21'
    containerOptions '--platform linux/amd64' 

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

    container 'rbceq2:2.4.1'

    publishDir {"results/${meta.id}"}, mode: 'copy' // after success, copy declared outputs into results/

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

    gvcf_ch = channel.fromPath(params.gvcfs)
        .map{gvcf -> [[id: gvcf.simpleName], gvcf]} // [[id:'HG00096'], file]
    
    CONVERT(gvcf_ch)

    RBCEQ2(CONVERT.out.vcf)
}