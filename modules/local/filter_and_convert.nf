process FILTER_AND_CONVERT {

    tag "${meta.id}"
     
    container params.bcftools_container

    input:
        tuple val(meta), path(gvcf)
    
    output:
        tuple val(meta), path("${meta.id}.converted.vcf.gz"), emit: vcf
    
    
    script:
    """
    bcftools norm -m -any -Ou $gvcf \
        | bcftools view -e 'ALT="<NON_REF>" || FMT/DP="." || FMT/DP < ${params.min_depth} || FMT/GQ="." || FMT/GQ < ${params.min_gq}' \
            --trim-alt-alleles -Oz -o ${meta.id}.converted.vcf.gz
    """
}