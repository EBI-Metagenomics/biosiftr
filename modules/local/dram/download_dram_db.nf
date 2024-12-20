process DOWNLOAD_DRAM_DB {

    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9'
        : 'biocontainers/gnu-wget:1.18--h36e9172_9'}"

    output:
    path("dram_dbs/"), emit: dram_db

    script:
    """
    mkdir -p dram_dbs

    wget -q --continue "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/v1.5.0/data/amg_database.tsv" -O dram_dbs/amg_database.tsv
    wget -q --continue "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/v1.5.0/data/etc_module_database.tsv" -O dram_dbs/etc_module_database.tsv
    wget -q --continue "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/v1.5.0/data/function_heatmap_form.tsv" -O dram_dbs/function_heatmap_form.tsv
    wget -q --continue "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/v1.5.0/data/genome_summary_form.tsv" -O dram_dbs/genome_summary_form.tsv
    wget -q --continue "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/v1.5.0/data/module_step_form.tsv" -O dram_dbs/module_step_form.tsv
    wget -q --continue "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz" -O dram_dbs/Pfam-A.hmm.dat.gz

    echo "Creating the CONFIG file for DRAM distill"
    cat > dram_dbs/DRAM_CONFIG.json << EOF
    {
        "description_db": "None",
        "kegg": null,
        "kofam": null,
        "kofam_ko_list": null,
        "uniref": null,
        "pfam": null,
        "pfam_hmm_dat": null,
        "dbcan": null,
        "dbcan_fam_activities": null,
        "viral": null,
        "peptidase": null,
        "vogdb": null,
        "vog_annotations": null,
        "genome_summary_form": "/data/genome_summary_form.tsv",
        "module_step_form": "/data/module_step_form.tsv",
        "etc_module_database": "/data/etc_module_database.tsv",
        "function_heatmap_form": "/data/function_heatmap_form.tsv",
        "amg_database": "/data/amg_database.tsv"
    }
    EOF
    """
}
