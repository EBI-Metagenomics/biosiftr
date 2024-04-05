#!/usr/bin/env python

import argparse
import os.path
import sys
import gzip
from Bio import SeqIO

##### This script transforms sourmash genomes relative abundance into species relative abundance
##### Alejandra Escobar, EMBL-EBI
##### Dec 20, 2023


def metadata_parser( catalogue_metadata ):
    ref_spec_genome = {}
    with open(catalogue_metadata, 'r') as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split('\t')
            rep_genome = l_line[13]
            lineage = l_line[14]+';'+rep_genome
            if not rep_genome in ref_spec_genome:
                ref_spec_genome[rep_genome] = lineage.replace(' ','_')

    return( ref_spec_genome )


def aggregate_species( sm_csv, ref_spec_genome, out_name ):
    prefix_name = out_name.replace('_species.tsv','')
    with gzip.open(sm_csv, 'rt', encoding='utf-8') as input_file, open(out_name,'w') as file_out:
        file_out.write("\t".join([
            'lineage',
            prefix_name
        ])+"\n")
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split(',')
            f_unique_weighted = l_line[4]
            name = l_line[9]

            if name in ref_spec_genome:
                lineage  = ref_spec_genome[name]
                file_out.write("\t".join([
                    lineage,
                    f_unique_weighted
                ])+"\n")
            else:
                print('No lineage for genome ',name)




def main():
    parser = argparse.ArgumentParser(
        description="This script transforms genomes relative abundance into species relative abundance. Provide the sourmash csv output file generated using kmer=51 vs the reps database and the genomes catalogue metadata file"
    )
    parser.add_argument(
        "--sm_csv",
        type=str,
        help="Sourmash output file in csv format (reps_gather-k51.csv)",
        required=True,
    )
    parser.add_argument(
        "--metadata",
        type=str,
        help="Catalogue metadata file",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Name of the output file. Default = sm_species_uw.tsv",
        required=False,
    )
    args = parser.parse_args()
    
    if args.output:
        out_name = args.output
    else:
        out_name = 'sm_species_uw.tsv'


    ### Calling functions
    ( ref_spec_genome ) = metadata_parser( args.metadata )
    aggregate_species( args.sm_csv, ref_spec_genome, out_name )

if __name__ == "__main__":
    main()

