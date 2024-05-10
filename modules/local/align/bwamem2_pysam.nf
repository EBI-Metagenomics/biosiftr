process ALIGN_BWAMEM2 {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/microbiome-informatics/bwa_eukcc:2.2.1_2.0'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.tsv"), emit: cov_file
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    echo " ---> mapping reads"
    bwa-mem2 \\
        mem \\
        -M \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | samtools view -@ ${task.cpus} -F256 -F4 -uS - \\
        | samtools sort -@ ${task.cpus} -O bam - -o ${prefix}_sorted.bam
    samtools index -@ ${task.cpus} ${prefix}_sorted.bam

    echo " ---> processing bam file"
    bam2cov_filt.py \\
        --bwa_bam ${prefix}_sorted.bam \\
        --prefix ${prefix}_u_relab_01

    echo " ---> removing bam file"
    rm -rf *.bam *.bai


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pysam: \$(python3 -c "import pkg_resources; print(pkg_resources.get_distribution('pysam').version)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_u_relab_01.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pysam: \$(python3 -c "import pkg_resources; print(pkg_resources.get_distribution('pysam').version)")
    END_VERSIONS
    """
}
