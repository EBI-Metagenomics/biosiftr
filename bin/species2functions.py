#!/usr/bin/env python

import argparse
import os.path
import sys
import gzip
from Bio import SeqIO

##### This script use the species prediction to generate functional tables from the pangenomic profiles
##### Alejandra Escobar, EMBL-EBI
##### Jan 11, 2024


def pfam_parser( pfam_data ):
    pfam_desc = {}
    with gzip.open(pfam_data, 'rt', encoding='utf-8') as input_file:
        entry = []
        for line in input_file:
            line = line.strip()
            if line == '//' and entry:
                entry_dict = {}
                for entry_line in entry:
                    if entry_line.startswith('#=GF DE'):
                        desc_string = entry_line.split('#=GF DE   ')[1]
                        entry_dict['DE'] = desc_string
                    elif entry_line.startswith('#=GF AC'):
                        entry_dict['AC'] = entry_line.split(' ')[-1].split('.')[0]
                if 'DE' in entry_dict and 'AC' in entry_dict:
                    pfam_desc[entry_dict['AC']] = entry_dict['DE']
                entry = []
            else:
                entry.append(line)
    return( pfam_desc )


def dram_parser( dram_form ):
    dram_desc = {}
    with open(dram_form, 'r') as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split('\t')
            gene_id = l_line[0]
            gene_description = l_line[1].replace('"','')
            if gene_id in dram_desc:
                if not gene_description in dram_desc[gene_id]:
                    dram_desc[gene_id].append(gene_description)
            else:
                dram_desc[gene_id] = [gene_description]
    return( dram_desc )


def relab_parser( relab_table ):
    taxonomy = {}
    reps_list = []
    with open(relab_table, 'r') as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split('\t')
            rep_genome = l_line[0].split(';')[-1]
            reps_list.append(rep_genome)
            lineage = l_line[0].split(';')
            lineage.pop(-1)
            lineage = [element for element in lineage if len(element) != 3]
            last_rank = lineage[-1]
            taxonomy[rep_genome] = last_rank
    return( reps_list, taxonomy )


def functions_finder_pan( reps_list, db_path ):
    functions_dict = {}
    species_kos = {}
    species_pfams = {}
    per_gene_dict = {}
    gene_positions = {}
    all_kos = []
    all_pfams = []
    for rep_genome in reps_list:
        db_file = db_path + '/' + rep_genome + '_clstr.tsv'
        positions = {}
        species_kos[rep_genome] = []
        species_pfams[rep_genome] = []
        per_gene_dict[rep_genome] = []
        pan_kos = []
        pan_pfams = []
        with open(db_file, 'r') as input_file:
            next(input_file)
            next(input_file)
            for line in input_file:
                per_gene_dict[rep_genome].append(line.rstrip())
                contig,gene_id,start,end,strand,kegg,cazy,pfam,core = line.rstrip().split('\t')

                if contig not in positions:
                    positions[contig] = {}
                    pos = 1
                else:
                    pos += 1
                gene_positions[gene_id] = str(pos)

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

    return( functions_dict, all_kos, all_pfams, species_kos, species_pfams, per_gene_dict, gene_positions )


def functions_finder_core( reps_list, db_path ):
    functions_dict = {}
    species_kos = {}
    species_pfams = {}
    per_gene_dict = {}
    gene_positions = {}
    all_kos = []
    all_pfams = []
    for rep_genome in reps_list:
        db_file = db_path + '/' + rep_genome + '_clstr.tsv'
        positions = {}
        species_kos[rep_genome] = []
        species_pfams[rep_genome] = []
        per_gene_dict[rep_genome] = []
        pan_kos = []
        pan_pfams = []
        with open(db_file, 'r') as input_file:
            next(input_file)
            next(input_file)
            for line in input_file:
                contig,gene_id,start,end,strand,kegg,cazy,pfam,core = line.rstrip().split('\t')

                if contig not in positions:
                    positions[contig] = {}
                    pos = 1
                else:
                    pos += 1
                gene_positions[gene_id] = str(pos)

                if core == 'true':
                    per_gene_dict[rep_genome].append(line.rstrip())
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

    return( functions_dict, all_kos, all_pfams, species_kos, species_pfams, per_gene_dict, gene_positions )


def community_writer( functions_dict, all_kos, all_pfams, output):
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


def dram_writer(per_gene_dict, gene_positions, taxonomy, pfam_desc, dram_desc, output):
    # The first column has the gene_id
    dram_header = [
        "",
        "fasta", 
        "scaffold", 
        "gene_position", 
        "start_position", 
        "end_position", 
        "strandedness", 
        "rank", 
        "kegg_id", 
        "kegg_hit", 
        "pfam_hits", 
        "cazy_hits", 
        "bin_taxonomy",
    ]

    with open(output+'_species_dram.tsv', 'w') as output_sp, \
        open(output+'_community_dram.tsv', 'w') as output_comm:
        output_sp.write('\t'.join(dram_header) + '\n')

        # Parsing the genes dictionary. We use the species_clstr id instead of fasta
        for species_clstr in per_gene_dict:
            for gene_line in per_gene_dict[species_clstr]:
                gene_info = {}

                # Populating handy info
                contig,gene_id,start,end,strand,kegg,pfam,cazy,core = gene_line.split('\t')
                print(core)
                gene_info["fasta"] = species_clstr+'_clstr: '+taxonomy[species_clstr]
                gene_info["scaffold"] = contig
                gene_info["start_position"] = start
                gene_info["end_position"] = end
                gene_info["strandedness"] = strand
                gene_info["bin_taxonomy"] = taxonomy[species_clstr]
                gene_info["gene_position"] = gene_positions[gene_id]

                # Processing cazy, pfam, and kegg descriptions. If no description,
                # then the annotation is discarded to avoid passing depricated annotation to DRAM
                rank = 'E'
                if cazy == '-':
                    cazy_hits = ''
                else:
                    cazy_desc_list = []
                    for cazy_acc in cazy.split(','):
                        if cazy_acc in dram_desc:
                            for cazy_desc in dram_desc[cazy_acc]:
                                if any([ cazy_desc.endswith(';'), cazy_desc.endswith('.') ]):
                                    cazy_desc = cazy_desc[:-1]
                                cazy_desc_list.append(cazy_desc)
                            last_element = cazy_desc_list.pop(-1)
                            last_element = last_element + ' [' + cazy_acc + ']'
                            cazy_desc_list.append(last_element)

                    if len(cazy_desc_list) > 0:
                        cazy_hits = '; '.join(cazy_desc_list)
                        rank = 'D'
                    else:
                        cazy_hits = ''
                gene_info["cazy_hits"] = cazy_hits

                if pfam == '-':
                    pfam_hits = ''
                else:
                    pfam_desc_list = []
                    for pfam_id in pfam.split(','):
                        if pfam_id in pfam_desc:
                            pfam_full_desc = pfam_desc[pfam_id]+' ['+pfam_id+']'
                            pfam_desc_list.append(pfam_full_desc)
                    if len(pfam_desc_list) > 0:
                        pfam_hits = ';'.join(pfam_desc_list)
                        rank = 'D'
                    else:
                        pfam_hits = ''
                gene_info["pfam_hits"] = pfam_hits

                if kegg == '-':
                    kegg_id = ''
                    kegg_hit = ''
                else:
                    kegg = kegg.replace(',',';')
                    kegg_id_list = []
                    ko_desc_list = []
                    for ko in kegg.split(';'):
                        if ko in dram_desc:
                            kegg_id_list.append(ko)
                            for ko_desc in dram_desc[ko]:
                                ko_desc_list.append(ko_desc)
                    if len(kegg_id_list) > 0:
                        rank = 'C'
                        kegg_id = ';'.join(kegg_id_list)
                        kegg_hit = ';'.join(ko_desc_list)
                    else:
                        kegg_id = ''
                        kegg_hit = ''
                gene_info["kegg_id"] = kegg_id
                gene_info["kegg_hit"] = kegg_hit
                gene_info["rank"] = rank

                # Writing to output files
                to_print = []
                to_print.append(gene_id)
                for header_key in dram_header[1:]:
                    to_print.append(gene_info[header_key])
                output_sp.write('\t'.join(to_print) + '\n')

                gene_info["fasta"] = 'community: '+output
                gene_info["bin_taxonomy"] = '' 
                to_print = []
                to_print.append(gene_id)
                for header_key in dram_header[1:]:
                    to_print.append(gene_info[header_key])
                output_comm.write('\t'.join(to_print) + '\n')
                output_sp.write('\t'.join(to_print) + '\n')



def main():
    parser = argparse.ArgumentParser(
        description="This script use the species prediction to generate functional tables from the pangenomic profiles"
    )
    parser.add_argument(
        "--pangenome_db",
        type=str,
        help="Path to the precomputed database of pangenomic functional profiles",
        required=True,
    )
    parser.add_argument(
        "--external_db",
        type=str,
        help="Path to the external db where Pfam-A.hmm.dat.gz AND genome_summary_form.tsv files exists",
        required=True,
    )
    parser.add_argument(
        "--relab",
        type=str,
        help="Species relative abundance table generated from BWA or Sourmash results. The first column is the representative species lineage and the corresponding genome ID",
        required=True,
    )
    parser.add_argument(
        "--core_mode",
        type=str,
        help="Either `core` or `pan` to indicate the fraction of the pangenome to be used for functional inference. The `core` mode is recommended for large catalogues like human-gut and mouse-gut",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Prefix to be used to name the output files. This string will also be used in headers",
        required=True,
    )
    args = parser.parse_args()
    
    pfam_db = args.external_db + "/Pfam-A.hmm.dat.gz"
    dram_form = args.external_db + "/genome_summary_form.tsv"

    ### Calling functions
    ( pfam_desc ) = pfam_parser( pfam_db )
    ( dram_desc ) = dram_parser( dram_form )
    ( reps_list, taxonomy ) = relab_parser( args.relab )

    if args.core_mode == 'pan':
        ( functions_dict, 
            all_kos, 
            all_pfams, 
            species_kos, 
            species_pfams,
            per_gene_dict,
            gene_positions ) = functions_finder_pan( reps_list, args.pangenome_db )
    elif args.core_mode == 'core':
        ( functions_dict,
            all_kos,
            all_pfams,
            species_kos,
            species_pfams,
            per_gene_dict,
            gene_positions ) = functions_finder_core( reps_list, args.pangenome_db )
    else:
        exit("The option "+args.core_mode+" is not allowed. Chose eaither `pan` or `core`")

    community_writer( 
            functions_dict, 
            all_kos, 
            all_pfams, 
            args.output)
    species_writer( 
            reps_list, 
            all_kos, 
            all_pfams, 
            species_kos, 
            species_pfams, 
            args.output)
    dram_writer( 
            per_gene_dict, 
            gene_positions, 
            taxonomy, 
            pfam_desc, 
            dram_desc, 
            args.output)


if __name__ == "__main__":
    main()

