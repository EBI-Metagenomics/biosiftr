#!/usr/bin/env python

import argparse
import os.path
import sys
from Bio import SeqIO

##### This script use the species prediction to generate functional tables from the pangenomic profiles
##### Alejandra Escobar, EMBL-EBI
##### Jan 11, 2024



def relab_parser( relab_table ):
    reps_list = []
    with open(relab_table, 'r') as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split('\t')
            rep_genome = l_line[0].split(';')[-1]
            reps_list.append(rep_genome)
    return( reps_list )


def functions_finder( reps_list, db_path ):
    functions_dict = {}
    species_kos = {}
    species_pfams = {}
    all_kos = []
    all_pfams = []
    for rep_genome in reps_list:
        db_file = db_path + '/' + rep_genome + '_clstr.tsv'
        species_kos[rep_genome] = []
        species_pfams[rep_genome] = []
        pan_kos = []
        pan_pfams = []
        with open(db_file, 'r') as input_file:
            next(input_file)
            next(input_file)
            for line in input_file:
                contig,gene_id,start,end,strand,kegg,pfam,cazy = line.rstrip().split('\t')
                if kegg != '-':
                    if ',' in kegg:
                        kegg_list = kegg.split(',')
                    else:
                        kegg_list = [kegg]
                    pan_kos = pan_kos + kegg_list
                    for current_ko in kegg_list:
                        if not current_ko in species_kos[rep_genome]:
                            species_kos[rep_genome].append(current_ko)

                if pfam != '-':
                    if ',' in pfam:
                        pfam_list = pfam.split(',')
                    else:
                        pfam_list = [pfam]
                    pan_pfams = pan_pfams + pfam_list
                    for current_pfam in pfam_list:
                        if not current_pfam in species_pfams[rep_genome]:
                            species_pfams[rep_genome].append(current_pfam)

        pan_kos = list(set(pan_kos))
        for ko in pan_kos:
            all_kos.append(ko)
            if ko in functions_dict:
                functions_dict[ko] += 1
            else:
                functions_dict[ko] = 1

        pan_pfams = list(set(pan_pfams))
        for pfam in pan_pfams:
            all_pfams.append(pfam)
            if pfam in functions_dict:
                functions_dict[pfam] += 1
            else:
                functions_dict[pfam] = 1

    all_kos = list(set(all_kos))
    all_pfams = list(set(all_pfams))

    return( functions_dict, all_kos, all_pfams, species_kos, species_pfams )


def community_writer( functions_dict, all_kos, all_pfams, output ):
    with open(output+'_community_kegg.tsv', 'w') as output_file:
        output_file.write('ko_id\t'+output+'\n')
        for ko in all_kos:
            output_file.write(ko+'\t'+str(functions_dict[ko])+'\n')

    with open(output+'_community_pfams.tsv', 'w') as output_file:
        output_file.write('pfam_id\t'+output+'\n')
        for pfam in all_pfams:
            output_file.write(pfam+'\t'+str(functions_dict[pfam])+'\n')


def species_writer( reps_list, all_kos, all_pfams, species_kos, species_pfams, output ):
    genomes_names = [genome + '_clstr' for genome in reps_list]

    with open(output+'_species_kegg.tsv', 'w') as output_file:
        header = ['ko_id'] + genomes_names
        output_file.write('\t'.join(header) + '\n')
        for ko in all_kos:
            to_print = []
            to_print.append(ko)
            for species in reps_list:
                if ko in species_kos[species]:
                    to_print.append('1')
                else:
                    to_print.append('0')
            output_file.write('\t'.join(to_print) + '\n')

    with open(output+'_species_pfams.tsv', 'w') as output_file:
        header = ['pfam_id'] + genomes_names
        output_file.write('\t'.join(header) + '\n')
        for pfam in all_pfams:
            to_print = []
            to_print.append(pfam)
            for species in reps_list:
                if pfam in species_pfams[species]:
                    to_print.append('1')
                else:
                    to_print.append('0')
            output_file.write('\t'.join(to_print) + '\n')



def main():
    parser = argparse.ArgumentParser(
        description="This script use the species prediction to generate functional tables from the pangenomic profiles"
    )
    parser.add_argument(
        "--db_path",
        type=str,
        help="Path to the precomputed database of pangenomic functional profiles",
        required=True,
    )
    parser.add_argument(
        "--relab",
        type=str,
        help="Species relative abundance table generated from BWA or Sourmash results. The first column is the representative species lineage and the corresponding genome ID",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Prefix to be used to name the output files. This string will also be used in headers",
        required=True,
    )
    args = parser.parse_args()
    
    ### Calling functions
    ( reps_list ) = relab_parser( args.relab )
    ( functions_dict, all_kos, all_pfams, species_kos, species_pfams ) = functions_finder( reps_list, args.db_path )
    community_writer( functions_dict, all_kos, all_pfams, args.output )
    species_writer( reps_list, all_kos, all_pfams, species_kos, species_pfams, args.output )

if __name__ == "__main__":
    main()

