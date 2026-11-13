process KEGG_PANGENOME_COMPLETENESS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kegg-pathways-completeness:1.4.4--pyhdfd78af_0':
        'biocontainers/kegg-pathways-completeness:1.4.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(profile)

    output:
    tuple val(meta), path("${meta.id}_clstr_kegg_comp.tsv"), emit: kegg_comp
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Column 6 of the profile is the KEGG KO list of a gene, column 9 the core flag.
    // The KO lists must have no empty entries and the tool must not run on an empty list:
    // kegg-pathways-completeness 1.4.4 counts an empty entry as a KO.
    """
    grep -v '^#' ${profile} | awk -F'\\t' '\$9 == "true" {print \$6}' \\
        | tr ',' '\\n' | grep -v -x -e '-' -e '' | sort -u | paste -sd, - > core_kos.txt
    grep -v '^#' ${profile} | awk -F'\\t' '{print \$6}' \\
        | tr ',' '\\n' | grep -v -x -e '-' -e '' | sort -u | paste -sd, - > pan_kos.txt

    tables=""
    for level in core pan; do
        if [ -s \${level}_kos.txt ]; then
            give_completeness -l \${level}_kos.txt -r \${level}
            tables="\${tables} --\${level} \${level}_pathways.tsv"
        fi
    done

    keggcomp_DB.py \${tables} --output ${meta.id}_clstr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg-pathways-completeness: \$(give_completeness --version | cut -d' ' -f3)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_clstr_kegg_comp.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kegg-pathways-completeness: \$(give_completeness --version | cut -d' ' -f3)
    END_VERSIONS
    """
}
