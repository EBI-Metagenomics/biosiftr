#!/bin/bash

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --biome)
            BIOME="$2"
            shift
            shift
            ;;
        --catalogue_dbs_path)
            CATALOGUE_DBS_PATH="$2"
            shift
            shift
            ;;
        --decont_refs_path)
            DECONT_REFS_PATH="$2"
            shift
            shift
            ;;
        --download_bwa)
            DOWNLOAD_BWA="$2"
            shift
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create verbose log file
LOG_FILE="dbs_setup_$(date +'%Y%m%d_%H%M%S').log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Change directory to decontamination references path
cd "$DECONT_REFS_PATH" || exit
if [ ! -d "reference_genomes" ]; then
    echo " ***  Creating the reference_genomes directory in $DECONT_REFS_PATH"
    mkdir reference_genomes && cd reference_genomes || exit
else
    echo " ***  The reference_genomes directory already exists in $DECONT_REFS_PATH"
    cd reference_genomes || exit
fi


# Check if human_phix.fa.* files exist
if ls human_phix.fa.* &>/dev/null; then
    echo " ***  The human and phiX reference genomes already exist. Skipping download"
else
    # Downloading human+phiX reference genomes
    echo " ***  Downloading the human and phiX reference genomes to ${DECONT_REFS_PATH}reference_genomes"
    wget --continue https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/human_phiX/human_phix_ref_bwamem2.tar.gz
    echo " ***  Extracting human and phiX reference genomes"
    tar -xvf human_phix_ref_bwamem2.tar.gz
    mv bwamem2/* .
    rm -r bwamem2 human_phix_ref_bwamem2.tar.gz
fi

# Check if $HOST.* files exist
HOST=$(echo "$BIOME" | cut -d '-' -f1)
if ls ${HOST}.* &>/dev/null; then
    echo " ***  The $HOST reference genome already exist. Skipping download"
else
    # Downloading the host genome
    echo " ***  Downloading the $HOST reference genome to $DECONT_REFS_PATH/reference_genomes"
    wget --continue "https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/$HOST/${HOST}_ref_bwamem2.tar.gz"
    echo " ***  Extracting the $HOST reference genome"
    tar -xvf "${HOST}_ref_bwamem2.tar.gz"
    mv bwamem2/* .
    rm -r bwamem2 "${HOST}_ref_bwamem2.tar.gz"
fi

# Downloading the catalogue-related files
cd "$CATALOGUE_DBS_PATH" || exit
if [ -d "$BIOME" ]; then
    echo " ***  A directory for the catalogue $BIOME already exists. Please remove the current directory to re-download. Exiting..."
    exit 1
else
    echo " ***  Creating $BIOME directory in $CATALOGUE_DBS_PATH"
    mkdir "$BIOME" && cd "$BIOME" || exit
fi

NEW_BIOME=$(echo $BIOME | sed 's/-vaginal-/-tmp-/;s/-v/|/;s/-tmp-/-vaginal-/' )
PREFIX_BIOME=$(echo "$NEW_BIOME" | cut -d '|' -f1)
VERSION=$(echo "$NEW_BIOME" | cut -d '|' -f2)
VERSION=$(echo "v$VERSION" | sed 's/-/./g' )

echo " ***  Downloading catalogue related databases to ${CATALOGUE_DBS_PATH}/${BIOME}"

# Downloading the catalogue metadata file
wget --continue "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$PREFIX_BIOME/$VERSION/genomes-all_metadata.tsv"


# Downloading the pangenome function tables
wget --continue "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$PREFIX_BIOME/$VERSION/pangenome_functions/functional_profiles.tar.gz"
tar -xvf functional_profiles.tar.gz
rm functional_profiles.tar.gz

wget --continue "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$PREFIX_BIOME/$VERSION/pangenome_functions/kegg_completeness.tar.gz"
tar -xvf kegg_completeness.tar.gz
rm kegg_completeness.tar.gz

# Downloading the representative genomes indexed for sourmash
wget --continue "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$PREFIX_BIOME/$VERSION/sourmash_db_${HOST}_${VERSION}/sourmash_species_representatives_k51.sbt.zip"

# Downloading bwamem2 db index if the option is set
if [ "$DOWNLOAD_BWA" = "true" ]; then
    echo " ***  Downloading bwamem2 indexed database for $BIOME to ${CATALOGUE_DBS_PATH}/${BIOME}"
    NEW_PREFIX=$(echo "$PREFIX_BIOME" | sed 's/-/_/')
    wget --continue "https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/${NEW_PREFIX}_reps/${NEW_PREFIX}-${VERSION}_bwamem2.tar.gz"
    tar -xvf "${NEW_PREFIX}-${VERSION}_bwamem2.tar.gz"
    mv "${NEW_PREFIX}-${VERSION}_bwamem2"/* .
    rm -r "${BIOME}-${VERSION}_bwamem2" "${NEW_PREFIX}-${VERSION}_bwamem2.tar.gz"
else
    echo " ***  Skipping download of bwamem2 indexed database for $BIOME"
    echo "      Note you will not be able to use --run_bwa true option on shallow-mapping pipeline for this biome"
fi

# Downloading external databases for dram visualization
cd "$CATALOGUE_DBS_PATH" || exit
if [ -d "external_dbs" ]; then
    echo " ***  Skipping external dbs downloading. The directory external_dbs already exists in $CATALOGUE_DBS_PATH"
else
    echo " ***  Downloading external dbs to $CATALOGUE_DBS_PATH/external_dbs/dram_distill_dbs"
    mkdir -p external_dbs/dram_distill_dbs && cd external_dbs/dram_distill_dbs || exit
    wget --continue "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/amg_database.tsv"
    wget --continue "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/etc_module_database.tsv"
    wget --continue "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/function_heatmap_form.tsv"
    wget --continue "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/genome_summary_form.tsv"
    wget --continue "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/module_step_form.tsv"
    wget --continue "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz"
fi

# Creating the CONFIG file for DRAM distill
echo " ***  Creating the CONFIG file for DRAM distill"
echo '{"description_db": "None", "kegg": null, "kofam": null, "kofam_ko_list": null, "uniref": null, "pfam": null, "pfam_hmm_dat": null, "dbcan": null, "dbcan_fam_activities": null, "viral": null, "peptidase": null, "vogdb": null, "vog_annotations": null, "genome_summary_form": "/data/genome_summary_form.tsv", "module_step_form": "/data/module_step_form.tsv", "etc_module_database": "/data/etc_module_database.tsv", "function_heatmap_form": "/data/function_heatmap_form.tsv", "amg_database": "/data/amg_database.tsv"}' > CONFIG


echo " ***  Databases setting up finished successfully for $BIOME"
echo " ***  Use the following parameters to run the shallow-mapping pipeline:"
echo "      nextflow run shallowmapping/main.nf \\"
echo "          --biome $BIOME \\"
echo "          --input your_samplesheet.csv \\"
echo "          --outdir your_outdir \\"
echo "          --shallow_dbs_path $CATALOGUE_DBS_PATH \\"
echo "          --decont_reference_paths ${DECONT_REFS_PATH}reference_genomes"




