process RBCEQ2 {

    tag "${meta.id}"

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
