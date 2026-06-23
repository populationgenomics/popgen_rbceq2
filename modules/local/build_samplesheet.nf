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