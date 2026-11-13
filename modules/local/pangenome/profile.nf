process PANGENOME_PROFILE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    // csv, prevalence, gff_dir and emapper are [] for singleton species
    tuple val(meta), path(rep_gff), path(rep_eggnog), path(csv), path(prevalence), path(gff_dir), path(emapper)
    path(pfam_dat)

    output:
    tuple val(meta), path("${meta.id}_clstr.tsv"), emit: profile
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def pangenome_args = csv ? "--csv ${csv} --prevalence ${prevalence} --gff_dir ${gff_dir}" : ''
    def emapper_arg = emapper ? "--emapper ${emapper}" : ''
    """
    pangenomeDB_builder_codon.py \\
        --mode post \\
        --rep ${meta.id} \\
        --rep_gff ${rep_gff} \\
        --rep_eggnog ${rep_eggnog} \\
        --pfam_dat ${pfam_dat} \\
        ${pangenome_args} \\
        ${emapper_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_clstr.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
