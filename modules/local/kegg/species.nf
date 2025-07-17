process KEGG_SPECIES {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(tax_table)
    val(tool)
    val(mode)
    path(kegg_db)

    output:
    tuple val(meta), path("*.tsv"), emit: spec_kegg_comp
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0' // WARN: Python script with no version control. This would be v1.0 of this script.
    """
    species2pathways.py \\
        --kegg_comp_db $kegg_db \\
        --relab $tax_table \\
        --core_mode $mode \\
        --output ${prefix}_${tool} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg_sp: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0' // WARN: Python script with no version control. This would be v1.0 of this script.
    """
    touch ${prefix}_${tool}_species_kegg_modules_comp.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg_sp: $VERSION
    END_VERSIONS
    """
}
