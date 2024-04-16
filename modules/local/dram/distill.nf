process DRAM_DISTILL {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dram:1.3.5--pyhdfd78af_0':
        'quay.io/biocontainers/dram:1.3.5--pyhdfd78af_0' }"

    containerOptions="--bind $params.shallow_dbs_path/external_dbs/dram_distill_dbs/:/data/ --bind $params.shallow_dbs_path/external_dbs/dram_distill_dbs/CONFIG:/usr/local/lib/python3.10/site-packages/mag_annotator/CONFIG"

    input:
    tuple val(meta), path(dram_summary)
    val(tool)       //sm for sourmash or bwa for bwamem2
    val(in_type)    //species or community

    output:
    tuple val(meta), path("*_dram*"), emit: destill_out
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.3.5' // WARN: dram has no option to print the tool version. This is the container version

    """
    if [[ "${in_type}" == "community" ]]; then
        echo ",fasta,scaffold,gene_position,start_position,end_position,strandedness,rank,kegg_id,kegg_hit,pfam_hits,cazy_hits,bin_taxonomy" | sed 's/,/\t/g' > dram_input.txt
    fi

    cat $dram_summary >> dram_input.txt  
    DRAM.py \\
        distill \\
        -i dram_input.txt  \\
        -o dram_out
    mv dram_out/product.html ${prefix}_${tool}_${in_type}_dram.html
    mv dram_out/product.tsv ${prefix}_${tool}_${in_type}_dram.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.3.5'
    """
    touch ${prefix}_${tool}_${in_type}_dram.html
    touch ${prefix}_${tool}_${in_type}_dram.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: $VERSION
    END_VERSIONS
    """
}
