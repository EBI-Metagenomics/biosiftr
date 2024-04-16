process KEGG_COMPLETENESS {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/microbiome-informatics/kegg-completeness:v1.1'

    input:
    tuple val(meta), path(kos_table)
    val(tool)
    
    output:
    tuple val(meta), path("*_community_kegg_modules_comp.tsv"), emit: kegg_comp
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1' // WARN: The tool have no option to print the running version. This is the container version
    """
    cat $kos_table | head -1 | sed 's/ko_id/module/' > ${prefix}_${tool}_community_kegg_modules_comp.tsv
    sed '1d' $kos_table | cut -f1 | tr "\n" ","  > kos.list

    run_pathways.sh \\
        -l kos.list \\
        -o result \\
        $args

    awk -F '\t' '{print \$1 "|" \$3 "\t" \$2}' result.summary.kegg_pathways.tsv | sed '/^module_accession/d' >> ${prefix}_${tool}_community_kegg_modules_comp.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg_comm: $VERSION
    END_VERSIONS

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1'
    """
    touch ${prefix}_${tool}_community_kegg_modules_comp.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg_comm: $VERSION
    END_VERSIONS
    """
}
