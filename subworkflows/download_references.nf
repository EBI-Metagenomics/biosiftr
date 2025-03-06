include { DOWNLOAD_DRAM_DB                      } from "../modules/local/dram/download_dram_db"
include { DOWNLOAD_HUMAN_PHIX_BWAMEM2_INDEX     } from '../modules/local/download_human_phix_bwamem2_index'
include { DOWNLOAD_MGNIFY_GENOMES_REFERENCE_DBS } from "../modules/local/download_mgnify_genomes_reference_dbs"

workflow DOWNLOAD_REFERENCES {
    take:
    biome
    bwamem2_mode
    host_name

    main:
    dram_dbs = Channel.empty()
    biome_sourmash_db = Channel.empty()
    biome_genomes_metadata = Channel.empty()
    biome_bwa_db = Channel.empty()
    biome_pangenome_functional_anns_db = Channel.empty()
    biome_kegg_completeness_db = Channel.empty()
    human_phix_index = Channel.empty()

    dram_dbs_dir = file("${params.reference_dbs}/dram_dbs")

    biome_folder = file("${params.reference_dbs}/${biome}")

    biome_sourmash_db_file = file("${params.reference_dbs}/${biome}/sourmash_species_representatives_k21.sbt.zip")
    biome_genomes_metadata_file = file("${params.reference_dbs}/${biome}/genomes-all_metadata.tsv")
    biome_bwa_db_files = file("${params.reference_dbs}/${biome}/bwa_reps.*")
    biome_pangenome_functional_anns_db_dir = file("${params.reference_dbs}/${biome}/functional_profiles_DB/")
    biome_kegg_completeness_db_dir = file("${params.reference_dbs}/${biome}/kegg_completeness_DB/")

    // DRAM databases
    dram_dbs = dram_dbs_dir.exists()
        ? dram_dbs_dir
        : DOWNLOAD_DRAM_DB().dram_db.first()

    // Human PhiX genome used to decontaminate the reads
    human_phix_index_dir = file("${params.decontamination_indexes}/human_phix*")
    human_phix_index = human_phix_index_dir.size() > 0
        ? human_phix_index_dir
        : DOWNLOAD_HUMAN_PHIX_BWAMEM2_INDEX().out.human_phix_index.first()

    // Check if all required files exist for the biome

    if ( biome_folder.exists() && biome_sourmash_db_file.exists() && biome_genomes_metadata_file.exists() && biome_pangenome_functional_anns_db_dir.exists() && biome_kegg_completeness_db_dir.exists() && (!bwamem2_mode || biome_bwa_db_files.size() > 0 )) {
        println "Genomes databases found in ${biome_folder}"
        biome_sourmash_db = biome_sourmash_db_file
        biome_genomes_metadata = biome_genomes_metadata_file
        biome_pangenome_functional_anns_db = biome_pangenome_functional_anns_db_dir
        biome_kegg_completeness_db = biome_kegg_completeness_db_dir
        biome_bwa_db = bwamem2_mode ? biome_bwa_db_files.collect() : Channel.empty()
    } else {
        println "Genomes databases NOT found. Running sbwf to download DBs"

        DOWNLOAD_MGNIFY_GENOMES_REFERENCE_DBS(biome, bwamem2_mode)
        biome_sourmash_db = DOWNLOAD_MGNIFY_GENOMES_REFERENCE_DBS.out.sourmash_db.first()
        biome_genomes_metadata = DOWNLOAD_MGNIFY_GENOMES_REFERENCE_DBS.out.genomes_metadata_tsv.first()
        biome_pangenome_functional_anns_db = DOWNLOAD_MGNIFY_GENOMES_REFERENCE_DBS.out.pangenome_functional_anns_db.first()
        biome_kegg_completeness_db = DOWNLOAD_MGNIFY_GENOMES_REFERENCE_DBS.out.kegg_completeness_db.first()
        biome_bwa_db = bwamem2_mode ? DOWNLOAD_MGNIFY_GENOMES_REFERENCE_DBS.out.bwamem2_index.collect() : Channel.empty()
    }

    human_phix_index_ch = channel.from( human_phix_index ).collect()
        .map { db_files ->
            [[id: "human_phix"], db_files]
        }

    emit:
    dram_dbs
    biome_sourmash_db
    biome_genomes_metadata
    biome_bwa_db
    biome_pangenome_functional_anns_db
    biome_kegg_completeness_db
    human_phix_index_ch
}
