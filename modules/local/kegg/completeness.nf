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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1' // WARN: The tool have no option to print the version. This is the container version
    """
    echo "creating kos.list file"
    cat ${kos_table} | awk 'NR>1' | cut -f1 | tr "\\n" "," > kos.list

    echo "Running KEGG completeness for the community"
    run_pathways.sh \\
        -l kos.list \\
        -o result

    echo "Formatting output"
    awk -v prefix=${prefix} \\
        -v tool=${tool} \\
        -F'\\t' \\
        'NR==1 {print "module\\t" prefix "_" tool; next} {print \$1 "|" \$3 "\\t" \$2}' \\
        result.summary.kegg_pathways.tsv > ${prefix}_${tool}_community_kegg_modules_comp.tsv

    echo "Reformatted output done"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg-pathways-completeness: $VERSION
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
        kegg-pathways-completeness: $VERSION
    END_VERSIONS
    """
}
