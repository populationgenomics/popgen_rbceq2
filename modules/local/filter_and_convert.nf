process FILTER_AND_CONVERT {

    tag "${meta.id}"
    
    container params.bcftools_container

    input:
        tuple val(meta), path(gvcf)
        val min_depth
        val min_gq
    
    output:
        tuple val(meta), path("${meta.id}.converted.vcf.gz"), emit: vcf
    
    /*
    We need to convert gVCFs into VCFs because the <NON_REF> ALT allele in gVCFs breaks RBCEQ2.
    Filter in bcftools, not in RBCeq2: RBCeq2 has no --depth or --quality parameters (its README was wrong)
    So depth/GQ filtering has to happen upstream here.

    1. bcftools norm -m -any: splits multiallelic records into one row per ALT allele.
    RBCeq2 looks at variants one allele at a time, so multiallelic sites have to be broken apart first or they get misread.

    2. bcftools view -e '...' — drops records matching any of these (||):
    - ALT="<NON_REF>": the gVCF <NON_REF> symbolic allele. Breaks RBCeq2.
    - FMT/DP="." || FMT/DP < min_depth — missing depth, or depth below threshold.
    - FMT/GQ="." || FMT/GQ < min_gq — missing genotype quality, or GQ below threshold.

    We explicitly test for "." and the < comparison because in bcftools a missing value doesn't reliably fail a < test.
    So without the "." check, no-data sites could slip through.

    3. --trim-alt-alleles: after filtering, removes ALT alleles no genotype references anymore.
    */
    script:
    """
    bcftools norm -m -any -Ou $gvcf \
        | bcftools view -e 'ALT="<NON_REF>" || FMT/DP="." || FMT/DP < ${min_depth} || FMT/GQ="." || FMT/GQ < ${min_gq}' \
            --trim-alt-alleles -Oz -o ${meta.id}.converted.vcf.gz
    """
}