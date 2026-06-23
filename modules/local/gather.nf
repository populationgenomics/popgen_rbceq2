process GATHER {
    
    tag "${cohort}"

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