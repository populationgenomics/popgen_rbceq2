process RBCEQ2 {

    tag "${meta.id}"

    container params.rbceq2_container

    input:
        tuple val(meta), path(vcf)
        val reference_genome
        val output_pdfs

    output:
        tuple val(meta), path("${meta.id}_*.tsv"), emit: results // globs all three .tsv files 
        path("${meta.id}_PDFs"), optional: true, emit: pdfs
    
    script:
    def make_pdfs = output_pdfs ? '--PDFs' : ''
    """
    rbceq2 \
        --vcf $vcf \
        --out ${meta.id} \
        --reference_genome ${reference_genome} \
        ${make_pdfs}
    """
}
