/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ebi-metagenomics/biosiftr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/ebi-metagenomics/biosiftr
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BIOSIFTR        } from './workflows/biosiftr'
include { BUILD_MGNIFY_DB } from './workflows/build_mgnify_db'

//
// WORKFLOW: Run main ebi-metagenomics/biosiftr analysis pipeline
//
workflow EBIMETAGENOMICS_BIOSIFTR {
    // Print help message if needed
    if (params.help) {
        def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
        def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
        def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --biome <BIOME> -profile docker\n" +
            "  or, to build the databases of a catalogue:\n" +
            "nextflow run ${workflow.manifest.name} --new_mgnify_db --catalogue_path <CATALOGUE_DIR> --eggnog_db_dir <EGGNOG_DIR> -profile docker"
        log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
        System.exit(0)
    }

    // Validate input parameters
    if (params.validate_params) {
        validateParameters()
    }

    WorkflowMain.initialise(workflow, params, log)

    if (params.new_mgnify_db) {
        BUILD_MGNIFY_DB ()
    } else {
        BIOSIFTR ()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    EBIMETAGENOMICS_BIOSIFTR ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
