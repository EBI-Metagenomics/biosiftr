process POSTPROC_BAM2COV {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.22.0--py39hcada746_0':
        'quay.io/biocontainers/pysam:0.22.0--py38h15b938a_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(cov)

    output:
    tuple val(meta), path("*.tsv"), emit: cov_file
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam2cov_filt.py \\
        --bwa_bam $bam \\
        --cov_thres $cov \\
        --prefix ${prefix}_u_relab \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pysam: \$(python3 -c "import pkg_resources; print(pkg_resources.get_distribution('pysam').version)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_u_relab.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pysam: \$(python3 -c "import pkg_resources; print(pkg_resources.get_distribution('pysam').version)")
    END_VERSIONS
    """
}
