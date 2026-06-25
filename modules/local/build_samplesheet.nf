process BUILD_SAMPLESHEET { 

    container params.metamist_container

    input:
        val project
        val cohorts

    output:
        path "*.tsv", emit: samplesheet

    script:
    def cohort_list = cohorts.join(' ')
    """
    fetch_cohort_samplesheet.py --project ${project} --cohorts ${cohort_list}
    """
    
}