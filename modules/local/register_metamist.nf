process REGISTER_METAMIST {

    tag "${cohort}"
    container params.metamist_container

    input:
        tuple val(cohort), path(tsvs)
        val project
        val outdir
        val output_version

    script:
    def out_prefix = "${outdir}/${output_version}/${cohort}/combined"
        """
        update_metamist.py \
            --project ${project} \
            --cohorts ${cohort} \
            --type custom \
            --output ${out_prefix}/combined.${cohort}.geno.tsv \
            --secondary pheno_numeric=${out_prefix}/combined.${cohort}.pheno_numeric.tsv pheno_alphanumeric=${out_prefix}/combined.${cohort}.pheno_alphanumeric.tsv

        """
}