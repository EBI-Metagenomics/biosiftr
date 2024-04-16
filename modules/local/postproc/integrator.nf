process POSTPROC_INTEGRATOR {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    path(files_list)
    val(annot_type)

    output:
    path("*.tsv")      , emit: integrated_matrix
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '1.0' // WARN: Python script with no version control. This would be v1.0 of this script.
    """
    matrix_integrator.py \\
        --input ${files_list.join(' ')} \\
        --output ${annot_type}_matrix.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        postproc: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def VERSION = '1.0'
    """
    touch ${annot_type}_matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        postproc: $VERSION
    END_VERSIONS
    """
}
