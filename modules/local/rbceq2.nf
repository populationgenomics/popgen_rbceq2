process RBCEQ2 {

    tag "${meta.id}"

    container params.rbceq2_container

    input:
        tuple val(meta), path(vcf)
        val reference_genome

    output:
        tuple val(meta), path("${meta.id}_*.tsv"), emit: results // globs all three .tsv files 
        path("${meta.id}_PDFs"), optional: true, emit: pdfs
    
    script:
    def args = task.ext.args ?: ''
    """
    rbceq2 \
        --vcf $vcf \
        --out ${meta.id} \
        --reference_genome ${reference_genome} \
        ${args}
    """
}
