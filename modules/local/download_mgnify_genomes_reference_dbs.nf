process DOWNLOAD_MGNIFY_GENOMES_REFERENCE_DBS {

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9'
        : 'biocontainers/gnu-wget:1.18--h36e9172_9'}"

    input:
    val biome
    val download_bwamem2

    output:
    tuple val(biome), path("genomes-all_metadata.tsv"),                     emit: genomes_metadata_tsv
    tuple val(biome), path("functional_profiles_DB/"),                      emit: pangenome_functional_anns_db
    tuple val(biome), path("kegg_completeness_DB/"),                        emit: kegg_completeness_db
    tuple val(biome), path("sourmash_species_representatives_k21.sbt.zip"), emit: sourmash_db
    tuple val(biome), path("bwamem2_index/")                              , emit: bwamem2_index, optional: true

    script:
    def matcher = biome =~ /(.+?)(-v[0-9.\.]+)?$/
    def biome_name = matcher[0][1]
    def biome_version = matcher[0][2] ? matcher[0][2].substring(1) : null
    if (!biome_version) {
        exit("Error the biome version of ${biome} can't be parsed.")
    }
    // MGnify genomes catalogue data //
    // Example: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/honeybee-gut/v1.0.1/
    def biome_catalogue_ftp = "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/${biome_name}/${biome_version}"

    // Shallow mapping specific //
    // This FTP path contains the MGnify Genomes catalogue processed annotations, ready to be used with this pipeline
    def ftp_base = "ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/mgnify_genomes/${biome_name}_reps/${biome_version}/"

    def functions_ftp = "${ftp_base}/pangenome_functional_profiles.tar.gz"
    def kegg_ftp = "${ftp_base}/kegg_completeness.tar.gz"
    def sourmash_ftp = "${ftp_base}/sourmash_species_representatives_k21.sbt.zip"
    def reps_bwamem2_index_ftp = "${ftp_base}/reps_bwamem2.tar.gz"

    """
    if [[ "${download_bwamem2}" == 'True' ]];
    then
        # Downloading the host genome #
        mkdir -p bwamem2_index/
        echo " ***  Downloading the biome reps mgnify genomes bwamem2 index ${reps_bwamem2_index_ftp}"
        wget -nv --show-progress --progress=bar:force:noscroll --continue "${reps_bwamem2_index_ftp}"

        echo " ***  Extracting the bwamem index..."
        tar -xvf reps_bwamem2.tar.gz -C bwamem2_index/
    fi

    # Downloading the catalogue-related files #
    echo " ***  Downloading catalogue related reference data"

    # Downloading the catalogue metadata file
    wget -nv --show-progress --progress=bar:force:noscroll --continue "${biome_catalogue_ftp}/genomes-all_metadata.tsv"

    wget -nv --show-progress --progress=bar:force:noscroll --continue "${functions_ftp}"
    tar -xvf pangenome_functional_profiles.tar.gz

    wget -nv --show-progress --progress=bar:force:noscroll --continue "${kegg_ftp}"
    tar -xvf kegg_completeness.tar.gz

    # Downloading the representative genomes indexed for sourmash
    wget -nv --show-progress --progress=bar:force:noscroll --continue "${sourmash_ftp}"
    """
}
