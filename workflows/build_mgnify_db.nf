/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BUILD_MGNIFY_DB: make the BioSIFTR databases of one MGnify genomes catalogue
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Input: --catalogue_path, the catalogue folder on the cluster
           (for example /nfs/production/rdf/metagenomics/catalogs/<host>/<vX.Y>/ftp)
    Output: <outdir>/<catalogue_id>/ with the layout that subworkflows/download_references.nf
            expects, the archives for the FTP site, and corrected pangenome files for the
            species that have duplicate gene names (catalogue_fixes/).
*/

include { paramsSummaryLog } from 'plugin/nf-validation'

include { PANGENOME_FIX                     } from '../modules/local/pangenome/fix'
include { CGT                               } from '../modules/local/pangenome/cgt'
include { PANGENOME_ACCESSORY               } from '../modules/local/pangenome/accessory'
include { PANGENOME_PROFILE                 } from '../modules/local/pangenome/profile'
include { KEGG_PANGENOME_COMPLETENESS       } from '../modules/local/kegg/pangenome_completeness'
include { DOWNLOAD_DRAM_DB                  } from '../modules/local/dram/download_dram_db'
include { EGGNOGMAPPER                      } from '../modules/nf-core/eggnogmapper/main'
include { SOURMASH_SKETCH as DB_SKETCH      } from '../modules/nf-core/sourmash/sketch/main'
include { SOURMASH_INDEX                    } from '../modules/nf-core/sourmash/index/main'
include { FIND_CONCATENATE as CAT_REPS      } from '../modules/nf-core/find/concatenate/main'
include { BWAMEM2_INDEX                     } from '../modules/nf-core/bwamem2/index/main'
include { TAR as TAR_FUNCTIONS              } from '../modules/nf-core/tar/main'
include { TAR as TAR_KEGG                   } from '../modules/nf-core/tar/main'
include { TAR as TAR_BWAMEM2                } from '../modules/nf-core/tar/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
    Catalogue ID from the catalogue path: .../<host>/v1.0/ftp -> <host>-v1-0
*/
def catalogueId(catalogue_dir) {
    def parts = catalogue_dir.toString().replaceAll('/+$', '').tokenize('/')
    if (parts[-1] == 'ftp') {
        parts = parts[0..-2]
    }
    def version = parts[-1]
    if (parts.size() < 2 || !(version ==~ /v[0-9]+(\.[0-9]+)*/)) {
        error("Cannot get the catalogue ID from '${catalogue_dir}'. Use --catalogue_id.")
    }
    return "${parts[-2]}-${version.replace('.', '-')}"
}

/*
    Folder prefix of a species: the accession without the version and the last 2 digits
    (MGYG000000761 -> MGYG0000007, MGYG000000761.1 -> MGYG0000007)
*/
def speciesPrefix(rep) {
    return rep.replaceAll(/\.[0-9]+$/, '')[0..-3]
}

/*
    A catalogue file can be compressed or not
*/
def catalogueFile(folder, name) {
    def plain = folder.resolve(name)
    def compressed = folder.resolve("${name}.gz")
    if (plain.exists()) {
        return plain
    }
    if (compressed.exists()) {
        return compressed
    }
    error("Catalogue file not found: ${plain}(.gz)")
}

workflow BUILD_MGNIFY_DB {

    log.info(paramsSummaryLog(workflow))

    ch_versions = channel.empty()

    catalogue_dir = file(params.catalogue_path, checkIfExists: true)
    metadata = file("${catalogue_dir}/genomes-all_metadata.tsv", checkIfExists: true)
    catalogue_id = params.catalogue_id ?: catalogueId(catalogue_dir)
    log.info("Catalogue ID: ${catalogue_id}")

    // The profiling workflow needs the catalogue metadata next to the databases
    channel.fromPath(metadata).collectFile(name: 'genomes-all_metadata.tsv', storeDir: "${params.outdir}/${catalogue_id}")

    eggnog_dmnd = file("${params.eggnog_db_dir}/eggnog_proteins.dmnd", checkIfExists: true)
    eggnog_data_dir = file(params.eggnog_db_dir, checkIfExists: true)

    // Pfam names to accessions, from the DRAM databases folder
    dram_dbs_dir = file("${params.reference_dbs}/dram_dbs")
    pfam_dat = dram_dbs_dir.exists()
        ? channel.value(file("${dram_dbs_dir}/Pfam-A.hmm.dat.gz", checkIfExists: true))
        : DOWNLOAD_DRAM_DB().dram_db.map { dram_db -> dram_db / 'Pfam-A.hmm.dat.gz' }.first()

    /*
        Species of the catalogue: [ meta, species folder ]
        The species folder is species_catalogue/<prefix>/<rep> (see speciesPrefix)
    */
    ch_species = channel.fromPath(metadata)
        .splitCsv(sep: '\t', header: true)
        .map { row -> [row.Species_rep, row.Genome] }
        .groupTuple()
        .map { rep, genomes ->
            def meta = [id: rep, prefix: speciesPrefix(rep), catalogue: catalogue_id, single: genomes.size() == 1]
            [meta, catalogue_dir.resolve("species_catalogue/${meta.prefix}/${rep}")]
        }
        .branch { _meta, _species_dir ->
            single: _meta.single
            pangenome: true
        }

    /*
        Stage A: check the pangenome files and correct the species with duplicate gene names
    */
    ch_pangenome_dirs = ch_species.pangenome.map { meta, species_dir ->
        [meta, species_dir.resolve('pan-genome'), catalogue_dir.resolve("all_genomes/${meta.prefix}/${meta.id}/genomes1")]
    }

    PANGENOME_FIX(ch_pangenome_dirs)
    ch_versions = ch_versions.mix(PANGENOME_FIX.out.versions.first())

    PANGENOME_FIX.out.report.collectFile(
        name: 'pangenome_check_report.tsv',
        keepHeader: true,
        skip: 1,
        sort: { report -> report.name },
        storeDir: "${params.outdir}/${catalogue_id}",
    )

    CGT(PANGENOME_FIX.out.fixed.join(PANGENOME_FIX.out.renamed), metadata)
    ch_versions = ch_versions.mix(CGT.out.versions.first())

    // Pangenome files of each species: [ meta, [csv, rtab, fasta, prevalence] ]
    // The corrected files replace the catalogue files when PANGENOME_FIX made them.
    ch_fixed_files = PANGENOME_FIX.out.fixed
        .join(CGT.out.cgt)
        .map { meta, csv, rtab, fasta, prevalence -> [meta, [csv, rtab, fasta, prevalence]] }

    ch_pangenome_files = ch_pangenome_dirs
        .map { meta, pangenome_dir, _gff_dir ->
            def files = [
                catalogueFile(pangenome_dir, 'gene_presence_absence.csv'),
                catalogueFile(pangenome_dir, 'gene_presence_absence.Rtab'),
                catalogueFile(pangenome_dir, 'pan-genome.fna'),
                catalogueFile(pangenome_dir, 'gene_prevalence_corrected.txt'),
            ]
            [meta, files]
        }
        .join(ch_fixed_files, remainder: true)
        .map { meta, catalogue_files, fixed_files -> [meta, fixed_files ?: catalogue_files] }

    /*
        Stage B: functional profiles
    */
    ch_gff_dirs = ch_pangenome_dirs.map { meta, _pangenome_dir, gff_dir -> [meta, gff_dir] }

    PANGENOME_ACCESSORY(ch_pangenome_files.join(ch_gff_dirs).map { meta, files, gff_dir -> [meta, files[0], gff_dir] })
    ch_versions = ch_versions.mix(PANGENOME_ACCESSORY.out.versions.first())

    EGGNOGMAPPER(
        PANGENOME_ACCESSORY.out.fasta.filter { _meta, fasta -> fasta.size() > 0 },
        ['diamond', eggnog_dmnd],
        eggnog_data_dir,
    )

    ch_rep_files = ch_species.single.mix(ch_species.pangenome).map { meta, species_dir ->
        def genome_dir = species_dir.resolve('genome')
        [meta, catalogueFile(genome_dir, "${meta.id}.gff"), catalogueFile(genome_dir, "${meta.id}_eggNOG.tsv")]
    }

    ch_profile_pangenome = ch_rep_files
        .filter { meta, _gff, _eggnog -> !meta.single }
        .join(ch_pangenome_files)
        .join(ch_gff_dirs)
        .join(EGGNOGMAPPER.out.annotations, remainder: true)
        .map { meta, gff, eggnog, files, gff_dir, emapper -> [meta, gff, eggnog, files[0], files[3], gff_dir, emapper ?: []] }

    ch_profile_single = ch_rep_files
        .filter { meta, _gff, _eggnog -> meta.single }
        .map { meta, gff, eggnog -> [meta, gff, eggnog, [], [], [], []] }

    PANGENOME_PROFILE(ch_profile_single.mix(ch_profile_pangenome), pfam_dat)
    ch_versions = ch_versions.mix(PANGENOME_PROFILE.out.versions.first())

    KEGG_PANGENOME_COMPLETENESS(PANGENOME_PROFILE.out.profile)
    ch_versions = ch_versions.mix(KEGG_PANGENOME_COMPLETENESS.out.versions.first())

    /*
        Stage B: indexes of the representative genomes
    */
    ch_rep_genomes = ch_species.single.mix(ch_species.pangenome).map { meta, species_dir ->
        [meta, catalogueFile(species_dir.resolve('genome'), "${meta.id}.fna")]
    }
    ch_db_meta = [id: catalogue_id, catalogue: catalogue_id]

    DB_SKETCH(ch_rep_genomes)
    ch_versions = ch_versions.mix(DB_SKETCH.out.versions.first())

    SOURMASH_INDEX(
        DB_SKETCH.out.signatures.map { _meta, signature -> signature }.collect(sort: true).map { signatures -> [ch_db_meta, signatures] },
        21,
    )

    CAT_REPS(
        ch_rep_genomes.map { _meta, genome -> genome }.collect(sort: true).map { genomes -> [ch_db_meta, genomes] }
    )

    BWAMEM2_INDEX(CAT_REPS.out.file_out, 'reps_bwamem2')

    /*
        Archives for the FTP site
    */
    TAR_FUNCTIONS(
        PANGENOME_PROFILE.out.profile.map { _meta, profile -> profile }.collect(sort: true).map { profiles -> [ch_db_meta + [id: 'pangenome_functional_profiles'], profiles] },
        '.gz',
    )
    TAR_KEGG(
        KEGG_PANGENOME_COMPLETENESS.out.kegg_comp.map { _meta, table -> table }.collect(sort: true).map { tables -> [ch_db_meta + [id: 'kegg_completeness'], tables] },
        '.gz',
    )
    TAR_BWAMEM2(
        BWAMEM2_INDEX.out.index.map { meta, index -> [meta + [id: 'reps_bwamem2'], index] },
        '.gz',
    )

    /*
        Versions (tools that use the 'versions' topic are added as versions.yml fragments)
    */
    ch_versions = ch_versions.mix(
        channel.topic('versions')
            .distinct()
            .map { process, tool, version -> "\"${process}\":\n  ${tool}: ${version}" }
            .collectFile(name: 'topic_versions.yml', newLine: true)
    )
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
}
