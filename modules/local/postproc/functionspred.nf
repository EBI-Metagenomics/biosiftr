process POSTPROC_FUNCTIONSPRED {
    tag "$meta.id"
    label 'process_single'


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(tax_tsv)
    val(tool)
    val(core_mode)
    path(pangenome_db)
    path(dram_dbs)
    val(run_dram)

    output:
    tuple val(meta), path("*_species_pfams.tsv")  , emit: pfam_spec
    tuple val(meta), path("*_community_pfams.tsv"), emit: pfam_comm
    tuple val(meta), path("*_species_kegg.tsv")   , emit: kegg_spec
    tuple val(meta), path("*_community_kegg.tsv") , emit: kegg_comm
    tuple val(meta), path("*_species_dram_summary.tsv")   , emit: dram_spec, optional: true
    tuple val(meta), path("*_community_dram_summary.tsv") , emit: dram_comm, optional: true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def dram_arg = (run_dram) ? "--dram_out" : ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0' // WARN: Python script with no version control. This would be v1.0 of this script.
    """
    species2functions.py \\
        --pangenome_db $pangenome_db \\
        --external_db $dram_dbs \\
        --relab $tax_tsv \\
        --core_mode $core_mode \\
        --output ${prefix}_${tool} \\
        ${dram_arg}

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
    touch ${tool}_${prefix}_community_dram_summary.tsv
    touch ${tool}_${prefix}_species_kegg.tsv
    touch ${tool}_${prefix}_species_pfams.tsv
    touch ${tool}_${prefix}_species_dram_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        postproc: $VERSION
    END_VERSIONS
    """
}
