process KEGG_COMPLETENESS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kegg-pathways-completeness:1.4.4--pyhdfd78af_0':
        'biocontainers/kegg-pathways-completeness:1.4.4--pyhdfd78af_0' }"

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
    """
    echo "creating kos.list file"
    # No trailing comma: kegg-pathways-completeness 1.4.4 counts an empty entry as a KO
    awk 'NR>1' ${kos_table} | cut -f1 | grep -v -x '' | paste -sd, - > kos.list

    echo "Running KEGG completeness for the community"
    if [ -s kos.list ]; then
        give_completeness \\
            -l kos.list \\
            -r result
    else
        # The tool reports a module for an empty list, so an empty list gives an empty table
        printf "module_accession\\tcompleteness\\tpathway_name\\n" > result_pathways.tsv
    fi

    echo "Formatting output"
    awk -v prefix=${prefix} \\
        -v tool=${tool} \\
        -F'\\t' \\
        'NR==1 {print "module\\t" prefix "_" tool; next} {print \$1 "|" \$3 "\\t" \$2}' \\
        result_pathways.tsv > ${prefix}_${tool}_community_kegg_modules_comp.tsv

    echo "Reformatted output done"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg-pathways-completeness: \$(give_completeness --version | cut -d' ' -f3)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${tool}_community_kegg_modules_comp.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg-pathways-completeness: \$(give_completeness --version | cut -d' ' -f3)
    END_VERSIONS
    """
}
