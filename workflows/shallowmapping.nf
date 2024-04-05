include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet; paramsSummaryMap } from 'plugin/nf-validation'

def summary_params = paramsSummaryMap(workflow)

validateParameters()

log.info paramsSummaryLog(workflow)

if (params.help) {
   log.info paramsHelp("nextflow run ebi-metagenomics-shallowmapping --help")
   exit 0
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


/*
~~~~~~~~~~~~~~~~~
    DBs
~~~~~~~~~~~~~~~~~
*/
dram_dbs          = file("$params.shallow_dbs_path/external_dbs/dram_distill_dbs", checkIfExists: true)
sourmash_db       = file("$params.shallow_dbs_path/$params.biome/sourmash_species_representatives_k51.sbt.zip", checkIfExists: true)
bwa_db            = file("$params.shallow_dbs_path/$params.biome/bwa_reps.fa.*", checkIfExists: true)
pangenome_db      = file("$params.shallow_dbs_path/$params.biome/functional_profiles_DB/", checkIfExists: true)
kegg_comp_db      = file("$params.shallow_dbs_path/$params.biome/kegg_completeness_DB/", checkIfExists: true)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Preprocessing modules
include { FASTP                               } from '../modules/local/fastp/main'
include { FASTQC as FASTQC_DECONT             } from '../modules/nf-core/fastqc/main'
include { MULTIQC                             } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS         } from '../modules/nf-core/custom/dumpsoftwareversions/main'

// Mapping modules
include { SOURMASH_GATHER                     } from '../modules/nf-core/sourmash/gather/main'
include { SOURMASH_SKETCH                     } from '../modules/nf-core/sourmash/sketch/main'
include { POSTPROC_SOURMASHTAXO               } from '../modules/local/postproc/sourmashtaxo'
include { POSTPROC_FUNCTIONSPRED as SM_FUNC   } from '../modules/local/postproc/functionspred'
include { DRAM_DISTILL as SM_DRAM             } from '../modules/local/dram/distill'
include { KEGG_COMPLETENESS as SM_COMM_KC     } from '../modules/local/kegg/completeness'
include { KEGG_SPECIES as SM_SPEC_KC          } from '../modules/local/kegg/species'

include { ALIGN_BWAMEM2                       } from '../modules/local/align/bwamem2'
include { POSTPROC_BAM2COV                    } from '../modules/local/postproc/bam2cov'
include { POSTPROC_BWATAXO                    } from '../modules/local/postproc/bwataxo'
include { POSTPROC_FUNCTIONSPRED as BWA_FUNC  } from '../modules/local/postproc/functionspred'
include { DRAM_DISTILL as BWA_DRAM            } from '../modules/local/dram/distill'
include { KEGG_COMPLETENESS as BWA_COMM_KC    } from '../modules/local/kegg/completeness'
include { KEGG_SPECIES as BWA_SPEC_KC         } from '../modules/local/kegg/species'

// Community results integration
include { POSTPROC_INTEGRATOR as INTEGRA_TAXO } from '../modules/local/postproc/integrator'
include { POSTPROC_INTEGRATOR as INTEGRA_KO   } from '../modules/local/postproc/integrator'
include { POSTPROC_INTEGRATOR as INTEGRA_PFAM } from '../modules/local/postproc/integrator'
include { POSTPROC_INTEGRATOR as INTEGRA_MODU } from '../modules/local/postproc/integrator'
include { DRAM_DISTILL as INTEGRA_DRAM        } from '../modules/local/dram/distill'

include { POSTPROC_INTEGRATOR as BWA_INT_TAXO } from '../modules/local/postproc/integrator'
include { POSTPROC_INTEGRATOR as BWA_INT_KO   } from '../modules/local/postproc/integrator'
include { POSTPROC_INTEGRATOR as BWA_INT_PFAM } from '../modules/local/postproc/integrator'
include { POSTPROC_INTEGRATOR as BWA_INT_MODU } from '../modules/local/postproc/integrator'
include { DRAM_DISTILL as BWA_INT_DRAM        } from '../modules/local/dram/distill'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READS_BWAMEM2_DECONTAMINATION as HUMAN_PHIX_DECONT } from '../subworkflows/ebi-metagenomics/reads_bwamem2_decontamination/main'
include { READS_BWAMEM2_DECONTAMINATION as HOST_DECONT       } from '../subworkflows/ebi-metagenomics/reads_bwamem2_decontamination/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def multiqc_report = []

workflow SHALLOWMAPPING {
    ch_versions       = Channel.empty()
    ch_log            = Channel.empty()
    ch_multiqc_files  = Channel.empty()


    // ---- Combine data into the reads channel ---- //
    groupReads = { list ->
        def meta = [id: list[0]]
        def reads = list.drop(1).findAll { it.size() > 0 }
        if (reads.size() == 1) {
            return tuple(meta + [single_end: true], reads)
        }
        else {
            return tuple(meta + [single_end: false], reads)
        }
    }
    ch_reads = Channel.fromSamplesheet("input").map(groupReads) // [ meta, [raw_reads] ]


    // ---- PREPROCESSING: Trimming, decontamination, and post-treatment qc ---- //
    FASTP ( ch_reads )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    // Creating channel for decontamination with human + phix genomes
    hp_ref = Channel.fromPath("$params.decont_reference_paths/human_phix*", checkIfExists: true).collect().map { db_files ->
        [ [id: 'human_phiX'], db_files ]
    }
    HUMAN_PHIX_DECONT ( FASTP.out.reads, hp_ref )
    ch_versions = ch_versions.mix(HUMAN_PHIX_DECONT.out.versions.first())

    // Creating channel for decontamination with host when biome != human
    def host_name = params.biome.split('-')[0]
    if ('human' in params.biome) {
        decont_reads = HUMAN_PHIX_DECONT.out.decontaminated_reads
    } else {
        host_ref = Channel.fromPath("$params.decont_reference_paths/$host_name.*", checkIfExists: true).collect().map { db_files ->
        [ [id: host_name], db_files ]
        }
        HOST_DECONT ( HUMAN_PHIX_DECONT.out.decontaminated_reads, host_ref )
        decont_reads = HOST_DECONT.out.decontaminated_reads
        ch_versions = ch_versions.mix(HOST_DECONT.out.versions.first())
    }


    // QC report after decontamination
    FASTQC_DECONT ( decont_reads )
    ch_versions = ch_versions.mix(FASTQC_DECONT.out.versions.first())


    // ---- MAPPING READS with sourmash: sketch decont reads, mapping, and profiling ---- //
    // Sketching decontaminated reads and running mapping
    SOURMASH_SKETCH ( decont_reads )
    ch_versions = ch_versions.mix(SOURMASH_SKETCH.out.versions.first())

    SOURMASH_GATHER ( SOURMASH_SKETCH.out.signatures, sourmash_db, false, false, false, false )
    ch_versions = ch_versions.mix(SOURMASH_GATHER.out.versions.first())

    // Processing sourmash mapping output: generating taxonomic and functional profiles
    POSTPROC_SOURMASHTAXO ( SOURMASH_GATHER.out.result, "$params.shallow_dbs_path/$params.biome/genomes-all_metadata.tsv" )
    ch_versions = ch_versions.mix(POSTPROC_SOURMASHTAXO.out.versions.first())    

    if (params.core_mode) {
        SM_FUNC ( POSTPROC_SOURMASHTAXO.out.sm_taxo, 'sm', 'core', pangenome_db, dram_dbs )
        ch_versions = ch_versions.mix(SM_FUNC.out.versions.first())

        SM_SPEC_KC ( POSTPROC_SOURMASHTAXO.out.sm_taxo, 'sm', 'core', kegg_comp_db )
        ch_versions = ch_versions.mix(SM_SPEC_KC.out.versions.first())
    } else {
        SM_FUNC ( POSTPROC_SOURMASHTAXO.out.sm_taxo, 'sm', 'pan', pangenome_db, dram_dbs )
        ch_versions = ch_versions.mix(SM_FUNC.out.versions.first())

        SM_SPEC_KC ( POSTPROC_SOURMASHTAXO.out.sm_taxo, 'sm', 'pan', kegg_comp_db )
        //ch_versions = ch_versions.mix(SM_SPEC_KC.out.versions.first())
    }

    SM_DRAM ( SM_FUNC.out.dram_spec, 'sm', 'species')
    ch_versions = ch_versions.mix(SM_DRAM.out.versions.first())

    SM_COMM_KC ( SM_FUNC.out.kegg_comm, 'sm' )
    //ch_versions = ch_versions.mix(SM_COMM_KC.out.versions.first())

    // ---- ANNOT INTEGRATOR: All samples matrices for taxo, kos, pfams, dram, and modules completeness ---- //
    INTEGRA_TAXO ( POSTPROC_SOURMASHTAXO.out.sm_taxo.collect{ it[1] }, 'sm_taxo' )
    ch_versions = ch_versions.mix(INTEGRA_TAXO.out.versions.first())

    INTEGRA_KO ( SM_FUNC.out.kegg_comm.collect{ it[1] }, 'sm_kos' )
    ch_versions = ch_versions.mix(INTEGRA_KO.out.versions.first())

    INTEGRA_PFAM ( SM_FUNC.out.pfam_comm.collect{ it[1] }, 'sm_pfam' )
    ch_versions = ch_versions.mix(INTEGRA_PFAM.out.versions.first())

    INTEGRA_MODU ( SM_COMM_KC.out.kegg_comp.collect{ it[1] }, 'sm_modules' )
    ch_versions = ch_versions.mix(INTEGRA_MODU.out.versions.first())

    //ch_dram_community = SM_FUNC.out.dram_comm.collectFile(name:'dram_community.tsv', newLine: true){ it[1] }.map { dram_summary -> [ [id: 'integrated'], dram_summary ] }
    //INTEGRA_DRAM ( ch_dram_community, 'sm', 'community' )
    //ch_versions = ch_versions.mix(INTEGRA_DRAM.out.versions.first())

    // ---- MAPPING READS with bwamem2 (optional): mapping, cleaning output, and profiling ---- //
    if (params.run_bwa) {
        genomes_ref = Channel.fromPath( bwa_db ).collect().map { db_files ->
        [ [id: host_name ], db_files ]
        }
        ALIGN_BWAMEM2 ( decont_reads, genomes_ref )
        ch_versions = ch_versions.mix(ALIGN_BWAMEM2.out.versions.first())

        POSTPROC_BAM2COV ( ALIGN_BWAMEM2.out.bam )
        ch_versions = ch_versions.mix(POSTPROC_BAM2COV.out.versions.first())

	POSTPROC_BWATAXO ( POSTPROC_BAM2COV.out.cov_file, "$params.shallow_dbs_path/$params.biome/genomes-all_metadata.tsv" )
	ch_versions = ch_versions.mix(POSTPROC_BWATAXO.out.versions.first())

        if (params.core_mode) {
            BWA_FUNC ( POSTPROC_BWATAXO.out.bwa_taxo, 'bwa', 'core', pangenome_db, dram_dbs )
            ch_versions = ch_versions.mix(BWA_FUNC.out.versions.first())

            BWA_SPEC_KC ( POSTPROC_BWATAXO.out.bwa_taxo, 'bwa', 'core', kegg_comp_db )
            //ch_versions = ch_versions.mix(BWA_SPEC_KC.out.versions.first())

        } else {
            BWA_FUNC ( POSTPROC_BWATAXO.out.bwa_taxo, 'bwa', 'pan', pangenome_db, dram_dbs )
            ch_versions = ch_versions.mix(BWA_FUNC.out.versions.first())

            BWA_SPEC_KC ( POSTPROC_BWATAXO.out.bwa_taxo, 'bwa', 'pan', kegg_comp_db )
            //ch_versions = ch_versions.mix(BWA_SPEC_KC.out.versions.first())
        }

        //BWA_DRAM (BWA_FUNC.out.dram_spec, 'bwa', 'species')
        //ch_versions = ch_versions.mix(BWA_DRAM.out.versions.first())

        BWA_COMM_KC ( BWA_FUNC.out.kegg_comm, 'bwa' )
        //ch_versions = ch_versions.mix(BWA_COMM_KC.out.versions.first())

        // ---- ANNOT INTEGRATOR: Matrices for taxo, kos, pfams, dram, and modules completeness ---- //
        BWA_INT_TAXO ( POSTPROC_BWATAXO.out.bwa_taxo.collect{ it[1] }, 'bwa_taxo' )
        ch_versions = ch_versions.mix(BWA_INT_TAXO.out.versions.first())

        BWA_INT_KO ( BWA_FUNC.out.kegg_comm.collect{ it[1] }, 'bwa_kos' )
        ch_versions = ch_versions.mix(BWA_INT_KO.out.versions.first())

        BWA_INT_PFAM ( BWA_FUNC.out.pfam_comm.collect{ it[1] }, 'bwa_pfam' )
        ch_versions = ch_versions.mix(BWA_INT_PFAM.out.versions.first())

        BWA_INT_MODU ( BWA_COMM_KC.out.kegg_comp.collect{ it[1] }, 'bwa_modules' )
        ch_versions = ch_versions.mix(BWA_INT_MODU.out.versions.first())

        ch_bwa_dram_community = BWA_FUNC.out.dram_comm.collectFile(name:'dram_community.tsv', newLine: true){ it[1] }.map { dram_summary -> [ [id: 'integrated'], dram_summary ] }
        //BWA_INT_DRAM ( ch_bwa_dram_community, 'bwa', 'community' )
        //ch_versions = ch_versions.mix(BWA_INT_DRAM.out.versions.first())
    }


    // ---- Multiqc report ---- //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    workflow_summary    = WorkflowShallowmapping.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowShallowmapping.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{ it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_DECONT.out.zip.collect{ it[1] }.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

}

