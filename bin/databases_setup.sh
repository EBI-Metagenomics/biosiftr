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
LOG_FILE="setup_log_$(date +'%Y%m%d_%H%M%S').log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Change directory to decontamination references path
cd "$DECONT_REFS_PATH" || exit
if [ ! -d "reference_genomes" ]; then
    mkdir reference_genomes && cd reference_genomes || exit
else
    cd reference_genomes || exit
fi

# Downloading human+phiX reference genomes
wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/human_phiX/human_phix_ref_bwamem2.tar.gz
tar -xvf human_phix_ref_bwamem2.tar.gz
mv bwamem2/* .
rm -r bwamem2

# Downloading the host genome
HOST=$(echo "$BIOME" | cut -d '-' -f1)
wget "https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/$HOST/${HOST}_ref_bwamem2.tar.gz"
tar -xvf "${HOST}_ref_bwamem2.tar.gz"
mv bwamem2/* .
rm -r bwamem2

# Downloading the catalogue-related files
cd "$CATALOGUE_DBS_PATH" || exit
if [ -d "$BIOME" ]; then
    echo "A database for the catalogue $BIOME already exists. Remove the current directory to re-download."
    exit 1
else
    mkdir "$BIOME" && cd "$BIOME" || exit
fi
VERSION=$(echo "$BIOME" | rev | cut -d '-' -f1 | rev)
VERSION="v$VERSION"
wget "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$HOST/$VERSION/genomes-all_metadata.tsv"

# Downloading the pangenome function tables
wget "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$HOST/$VERSION/pangenome_functions/functional_profiles.tar.gz"
wget "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$HOST/$VERSION/pangenome_functions/kegg_completeness.tar.gz"
tar -xvf kegg_completeness.tar.gz

# Downloading the representative genomes indexed for sourmash
wget "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$HOST/$VERSION/sourmash_db_${HOST}_${VERSION}/sourmash_species_representatives_k51.sbt.zip"

# Downloading bwamem2 db index if the option is set
if [ "$DOWNLOAD_BWA" = "true" ]; then
    wget "https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/${BIOME}_reps/${BIOME}-${VERSION}_bwamem2.tar.gz"
    tar -xvf "${BIOME}-${VERSION}_bwamem2.tar.gz"
    mv "${BIOME}-${VERSION}_bwamem2"/* .
    rm -r "${BIOME}-${VERSION}_bwamem2"
fi

# Downloading external databases for dram visualization
cd "$CATALOGUE_DBS_PATH" || exit
if [ -d "external_dbs" ]; then
    echo "Skipping external dbs downloading. The directory external_dbs already exists in $CATALOGUE_DBS_PATH."
else
    mkdir -p external_dbs/dram_distill_dbs && cd external_dbs/dram_distill_dbs || exit
    wget "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/amg_database.tsv"
    wget "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/etc_module_database.tsv"
    wget "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/function_heatmap_form.tsv"
    wget "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/genome_summary_form.tsv"
    wget "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/module_step_form.tsv"
    wget "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz"
fi

echo "Databases setting up finished successfully for $BIOME"

