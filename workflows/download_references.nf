include { DOWNLOAD_REFERENCES } from '../subworkflows/download_references'

workflow {
    DOWNLOAD_REFERENCES(params.biome, params.run_bwa)
}
