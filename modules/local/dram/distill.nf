process DRAM_DISTILL {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dram:1.3.5--pyhdfd78af_0':
        'quay.io/biocontainers/dram:1.3.5--pyhdfd78af_0' }"

    containerOptions {
        def arg = ""
        switch (workflow.containerEngine) {
            case 'singularity':
                arg = "--bind"
                break;
            case 'docker':
                arg = "--volume"
                break;
        }
        mounts = [
            "${params.shallow_dbs_path}/external_dbs/dram_distill_dbs/:/data/",
            "${params.shallow_dbs_path}/external_dbs/dram_distill_dbs/CONFIG:/usr/local/lib/python3.10/site-packages/mag_annotator/CONFIG"
        ]
        return "${arg} " + mounts.join(" ${arg} ")
    }

    input:
    tuple val(meta), path(dram_summary)
    val(tool)       //sm for sourmash or bwa for bwamem2
    val(in_type)    //species or community

    output:
    tuple val(meta), path("*_dram*"), emit: destill_out, optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.3.5' // WARN: dram has no option to print the tool version. This is the container version

    """
    if [[ "${in_type}" == "community" ]]; then
        echo ",fasta,scaffold,gene_position,start_position,end_position,strandedness,rank,kegg_id,kegg_hit,pfam_hits,cazy_hits,bin_taxonomy" | sed 's/,/\t/g' > dram_input.tsv
    fi

    cat $dram_summary >> dram_input.tsv
    line_count=\$(wc -l dram_input.tsv | cut -d' ' -f1)

    echo "Line count is "\$line_count

    # Just in case, remove the folder
    rm -rf dram_out/

    if [[ \$line_count > 1 ]]; then
        DRAM.py \\
            distill \\
            -i dram_input.tsv  \\
            -o dram_out
        mv dram_out/product.html ${prefix}_${tool}_${in_type}_dram.html
        mv dram_out/product.tsv ${prefix}_${tool}_${in_type}_dram.tsv
    else
        echo "The dram_input.tsv file is empty... skpping dram"
    fi

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
