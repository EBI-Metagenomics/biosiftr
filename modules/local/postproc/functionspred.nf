process POSTPROC_FUNCTIONSPRED {
    tag "$meta.id"
    label 'process_low'


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'quay.io/biocontainers/biopython:1.81' }"


    input:
    tuple val(meta), path(tax_tsv)
    path(pangenome_db)
    val(tool)

    output:
    tuple val(meta), path("*.tsv"), emit: func_tables
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0' // WARN: Python script with no version control. This would be v1.0 of this script.
    """
    species2functions.py \\
        --db_path $pangenome_db \\
        --relab $tax_tsv \\
        --output ${prefix}_${tool} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        postproc: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${tool}_${prefix}_community_kegg.tsv
    touch ${tool}_${prefix}_community_pfams.tsv
    touch ${tool}_${prefix}_species_kegg.tsv
    touch ${tool}_${prefix}_species_pfams.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        postproc: $VERSION
    END_VERSIONS
    """
}
