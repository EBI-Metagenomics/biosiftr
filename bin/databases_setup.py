#!/usr/bin/env python

import argparse
import os.path
import sys

##### This script download the databases for the shallow-mapping pipeline
##### Alejandra Escobar, EMBL-EBI
##### April 3, 2024

def core_parser( core_table ):
    all_modules = []
    core_values = {}
    with open( core_table, 'r' ) as input_file:
        next( input_file )
        for line in input_file:
            l_line = line.rstrip().split('\t')
            module_accession = l_line[0]
            completeness = l_line[1]
            pathway_name = l_line[2]
            module_name = module_accession + '|' + pathway_name
            all_modules.append(module_name)
            core_values[module_name]=completeness
    return( all_modules, core_values )


def pan_parser( all_modules, pan_table ):
    pan_values = {}
    with open( pan_table, 'r' ) as input_file:
        next( input_file )
        for line in input_file:
            l_line = line.rstrip().split('\t')
            module_accession = l_line[0]
            completeness = l_line[1]
            pathway_name = l_line[2]
            module_name = module_accession + '|' + pathway_name
            all_modules.append(module_name)
            pan_values[module_name]=completeness
    all_modules = list(set(all_modules))
    return( all_modules, pan_values )


def output_writer( all_modules, core_values, pan_values, prefix ):
    with open(prefix+'_kegg_comp.tsv','w') as file_out:
        file_out.write("\t".join([
            '#module',
            'pangenome',
            'core',
        ])+'\n')
        for module in all_modules:
            to_print = []
            to_print.append(module)
            if module in pan_values:
                to_print.append(pan_values[module])
            else:
                to_print.append('0')
            
            if module in core_values:
                to_print.append(core_values[module])
            else:
                to_print.append('0')

            file_out.write("\t".join(to_print)+'\n')



def main():
    parser = argparse.ArgumentParser(
        description="This script download the databases for the shallow-mapping pipeline"
    )
    parser.add_argument(
        "--biome",
        type=str,
        help="Provide the catalogue ID name",
        required=True,
    )
    parser.add_argument(
        "--catalogue_dbs_path",
        type=str,
        help="Location of the catalogue databases",
        required=True,
    )
    parser.add_argument(
        "--decont_refs_path",
        type=str,
        help="Location of the decontamination reference genomes",
        required=True,
    )
    parser.add_argument(
        "--download_bwa",
        type=bool,
        help="Download bwamem2 indexed reference genomes <default = false>",
        required=True,
    )
    args = parser.parse_args()
    
    ### Calling functions
    ( all_modules, core_values ) = core_parser( args.core )
    ( all_modules, pan_values ) = pan_parser( all_modules, args.pan )
    output_writer( all_modules, core_values, pan_values, args.output )

'''
mkdir -p databases/reference_genomes && cd databases/reference_genomes

# Downloading human+phiX reference genomes

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/human_phiX/human_phix_ref_bwamem2.tar.gz
tar -xvf human_phix_ref_bwamem2.tar.gz
mv bwamem2/* .
rm -r bwamem2

# Downloading the host genome. Replace $HOST by the name of the reference genome you are intend to download (e.g. chocken)

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/$HOST/$HOST_ref_bwamem2.tar.gz
tar -xvf $HOST_ref_bwamem2.tar.gz
mv bwamem2/* .
rm -r bwamem2

cd shallowmapping/databases/
mkdir $CATALOGUE_ID && cd $CATALOGUE_ID

# Downloading the catalogue metadata file. Replace $HOST for the name of the catalogue and $VERSION for the version you are intend to use as a reference

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$HOST/$VERSION/genomes-all_metadata.tsv 

# Downloading the pangenome function tables

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$HOST/$VERSION/pangenome_functions/functional_profiles.tar.gz
tar -xvf functional_profiles.tar.gz

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$HOST/$VERSION/pangenome_functions/kegg_completeness.tar.gz
tar -xvf kegg_completeness.tar.gz

# Downloading the representative genomes indexed for sourmash

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/$HOST/$VERSION/sourmash_db_$HOST_$VERSION/sourmash_species_representatives_k51.sbt.zip


wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/references/$BIOME_reps/$BIOME-$VERSION_bwamem2.tar.gz
tar -cvf $BIOME-$VERSION_bwamem2.tar.gz
mv $BIOME-$VERSION_bwamem2/* .
rm -r $BIOME-$VERSION_bwamem2


cd shallowmapping/databases/
mkdir -p external_dbs/dram_distill_dbs && cd external_dbs/dram_distill_dbs
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/amg_database.tsv
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/etc_module_database.tsv
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/function_heatmap_form.tsv
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/genome_summary_form.tsv
wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM/master/data/module_step_form.tsv
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz




'''



if __name__ == "__main__":
    main()

