process PANGENOME_FIX {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(pangenome_dir), path(gff_dir)

    output:
    tuple val(meta), path("fixed/gene_presence_absence.csv.gz"), path("fixed/gene_presence_absence.Rtab.gz"), path("fixed/pan-genome.fna.gz"), emit: fixed, optional: true
    tuple val(meta), path("fixed/renamed_genes.tsv")                                                                                        , emit: renamed, optional: true
    path "${meta.id}_pangenome_check.tsv"                                                                                                   , emit: report
    path "versions.yml"                                                                                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pangenome_fix.py \\
        --pangenome_dir ${pangenome_dir} \\
        --gff_dir ${gff_dir} \\
        --species ${meta.id} \\
        --outdir fixed \\
        --report ${meta.id}_pangenome_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    printf "species\\tstatus\\n${meta.id}\\tok\\n" > ${meta.id}_pangenome_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
