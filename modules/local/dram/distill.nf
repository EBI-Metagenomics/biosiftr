process DRAM_DISTILL {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dram:1.4.6--pyhdfd78af_2':
        'quay.io/biocontainers/dram:1.4.6--pyhdfd78af_2' }"

    containerOptions="--bind $params.dram_dbs/:/data/ --bind $params.dram_dbs/CONFIG:/usr/local/lib/python3.11/site-packages/mag_annotator/CONFIG"


    input:
    tuple val(meta), path(dram_summary)
    val(tool)       //sm for sourmash or bwa for bwamem2
    val(in_type)    //species or community

    output:
    tuple val(meta), path("dram_out/*"), emit: destill_out
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.4.6' // WARN: dram has no option to print the tool version. This is the container version
    """
    DRAM.py \\
        distill \\
        -i $dram_summary  \\
        -o dram_out
    mv dram_out/genome_stats.tsv dram_out/${prefix}_${tool}_${in_type}_stats.tsv
    mv dram_out/metabolism_summary.xlsx dram_out/${prefix}_${tool}_${in_type}_summary.xlsx
    mv dram_out/product.html dram_out/${prefix}_${tool}_${in_type}.html
    mv dram_out/product.tsv dram_out/${prefix}_${tool}_${in_type}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.4.6' // WARN: dram has no option to print the version of the tool
    """
    mkdir dram_out
    touch dram_out/${prefix}_${tool}_${in_type}_stats.tsv
    touch dram_out/${prefix}_${tool}_${in_type}_summary.xlsx
    touch dram_out/${prefix}_${tool}_${in_type}.html
    touch dram_out/${prefix}_${tool}_${in_type}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: $VERSION
    END_VERSIONS
    """
}
