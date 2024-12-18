process DOWNLOAD_HUMAN_PHIX_BWAMEM2_INDEX {

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9'
        : 'biocontainers/gnu-wget:1.18--h36e9172_9'}"

    output:
    path ("human_phix.fa*"), emit: human_phix_index

    script:
    """
    wget -q --continue ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/human_phiX/human_phix_ref_bwamem2.tar.gz
    tar -xvf human_phix_ref_bwamem2.tar.gz
    """
}
