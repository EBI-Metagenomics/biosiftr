process PANGENOME_ACCESSORY {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(csv), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_accessory.fasta"), emit: fasta
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pangenomeDB_builder_codon.py \\
        --mode pre \\
        --rep ${meta.id} \\
        --csv ${csv} \\
        --fasta ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_accessory.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
