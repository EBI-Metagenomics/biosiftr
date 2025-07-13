process DRAM_DISTILL {
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/dram:1.3.5--pyhdfd78af_0'
        : 'quay.io/biocontainers/dram:1.3.5--pyhdfd78af_0'}"

    containerOptions {
        def arg = "--volume"
        if (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') {
            arg = "--bind"
        }
        def mounts = [
            "${task.workDir}/${reference_dbs}/dram_dbs/:/data/",
            "${task.workDir}/${reference_dbs}/dram_dbs/DRAM_CONFIG.json:/usr/local/lib/python3.10/site-packages/mag_annotator/CONFIG",
        ]
        return "${arg} " + mounts.join(" ${arg} ")
    }

    input:
    tuple val(meta), path(dram_summary)
    path(reference_dbs)
    val tool      // sm for sourmash or bwa for bwa-mem2
    val in_type   // species or community

    output:
    tuple val(meta), path("*_dram*"), emit: destill_out, optional: true
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // WARN: dram has no option to print the tool version. This is the container version
    def VERSION = '1.3.5'
    def header_columns = [
        "",
        "fasta",
        "scaffold",
        "gene_position",
        "start_position",
        "end_position",
        "strandedness",
        "rank",
        "kegg_id",
        "kegg_hit",
        "pfam_hits",
        "cazy_hits",
        "bin_taxonomy"
    ]
    header_columns = header_columns.join("\t")
    """
    rm -rf dram_out/

    if [[ "${in_type}" == "community" ]]; then
        echo "${header_columns}" > dram_input.tsv
    fi

    cat ${dram_summary} >> dram_input.tsv

    line_count=\$(wc -l dram_input.tsv | cut -d' ' -f1)

    echo "Line count is "\$line_count

    if [[ \$line_count > 1 ]]; then
        DRAM.py \\
            distill \\
            -i dram_input.tsv  \\
            -o dram_out

        export counter=0
        # Loop through each product_*.html files #
        for productfile in dram_out/product*.html; do
            mv "\$productfile" "${prefix}_${tool}_${in_type}_\${counter}_dram.html"
            counter=\$((counter + 1))
        done

        if [ -f dram_out/product.html ]; then
            mv "dram_out/product.html ${prefix}_${tool}_${in_type}_dram.html"
        fi

        mv dram_out/product.tsv ${prefix}_${tool}_${in_type}_dram.tsv
    else
        echo "The dram_input.tsv file is empty... skpping dram"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: ${VERSION}
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
        dram: ${VERSION}
    END_VERSIONS
    """
}
