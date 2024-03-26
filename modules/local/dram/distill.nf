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
    tuple val(meta), path("*_dram.{html,tsv}"), emit: destill_out
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.4.6' // WARN: dram has no option to print the tool version. This is the container version

    """
    if [[ "${in_type}" == "community" ]]; then
        echo ",fasta,scaffold,gene_position,start_position,end_position,strandedness,rank,kegg_id,kegg_hit,pfam_hits,cazy_best_hit,bin_taxonomy" | sed 's/,/\t/g' > community_input.txt
        cat $dram_summary >> community_input.txt  
        DRAM.py \\
            distill \\
            -i community_input.txt  \\
            -o dram_out
        mv dram_out/product.html ${prefix}_${tool}_${in_type}_dram.html
        mv dram_out/product.tsv ${prefix}_${tool}_${in_type}_dram.tsv
    else
        DRAM.py \\
            distill \\
            -i $dram_summary  \\
            -o dram_out
        mv dram_out/product.html ${prefix}_${tool}_${in_type}_dram.html
        mv dram_out/product.tsv ${prefix}_${tool}_${in_type}_dram.tsv
    fi

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
    touch ${prefix}_${tool}_${in_type}_stats.tsv
    touch ${prefix}_${tool}_${in_type}_summary.xlsx
    touch ${prefix}_${tool}_${in_type}.html
    touch ${prefix}_${tool}_${in_type}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: $VERSION
    END_VERSIONS
    """
}
