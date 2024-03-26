#!/usr/bin/env python

import argparse
import os.path
import sys
import wget
import gzip
from Bio import SeqIO

##### This script find the accessory genes that needs eggNOG annotation
##### Alejandra Escobar, EMBL-EBI
##### Jan 8, 2024


def pfam_parser( pfam_data ):
    pfam_desc = {}
    with gzip.open(pfam_data, 'rt', encoding='utf-8') as input_file:
        entry = []
        for line in input_file:
            line = line.strip()
            if line == '//' and entry:
                entry_dict = {}
                for entry_line in entry:
                    if entry_line.startswith('#=GF ID'):
                        entry_dict['ID'] = entry_line.split()[-1]
                    elif entry_line.startswith('#=GF AC'):
                        entry_dict['AC'] = entry_line.split()[-1].split('.')[0]
                if 'ID' in entry_dict and 'AC' in entry_dict:
                    pfam_desc[entry_dict['ID']] = entry_dict['AC']
                entry = []
            else:
                entry.append(line)
    return( pfam_desc )


def metadata_parser( catalogue_metadata ):
    reps_clusters = {}
    with open( catalogue_metadata, 'r' ) as input_file:
        next( input_file )
        for line in input_file:
            l_line = line.rstrip().split('\t')
            clstr_member = l_line[0]
            rep_genome = l_line[13]
            if rep_genome in reps_clusters:
                reps_clusters[rep_genome].append(clstr_member)
            else:
                reps_clusters[rep_genome] = [clstr_member]
    return( reps_clusters )


def accessory_writer( reps_clusters, url_prefix ):
    for rep in reps_clusters:
        if len(reps_clusters[rep]) > 1:
            rep_prefix = rep[:-2]
            pan_url = url_prefix+'/species_catalogue/'+rep_prefix+'/'+rep+'/pan-genome/'

            # Downloading and parsing the presence/absence tab
            r_tab_url = pan_url+'gene_presence_absence.Rtab'
            r_tab_out = rep+'_pan.Rtab'
            if not os.path.exists(r_tab_out):
                try:
                    wget.download(r_tab_url, out=r_tab_out)
                except Exception as e:
                    print(f"Failed to download the file: {str(e)}")

            accessory_genes = []
            with open(r_tab_out, 'r') as input_file:
                header = input_file.readline().strip().split('\t')
                index = header.index(rep)
                for line in input_file:
                    l_line = line.rstrip().split('\t')
                    if l_line[index] == '0':
                        accessory_genes.append(l_line[0])
            os.remove(r_tab_out)

            # Downloading and parsing the fasta file of genes
            genes_url = pan_url+'pan-genome.fna'
            genes_out = rep+'_pan.fna'
            if not os.path.exists(genes_out):
                try:
                    wget.download(genes_url, out=genes_out)
                except Exception as g:
                    print(f"Failed to download the file: {str(g)}")

            with open(rep+'_accessory.fasta', 'w') as fasta_out:
                for record in SeqIO.parse(genes_out, "fasta"):
                    seq_id = str(record.id)
                    if seq_id in accessory_genes:
                        fasta_out.write('>'+seq_id+'\n')
                        fasta_out.write(str(record.seq).upper()+'\n')
            os.remove(genes_out)


def annot_writer( reps_clusters, url_prefix, pfam_desc ):
    for rep in reps_clusters:
        core_list, core_mgygs = [], []
        rep_prefix = rep[:-2]
        rep_url = url_prefix+'/species_catalogue/'+rep_prefix+'/'+rep+'/genome/'

        # Downloading and parsing the gff file of clusters size =1
        gff_url = rep_url+rep+'.gff'
        gff_out = rep+'.gff'
        if not os.path.exists(gff_out):
            try:
                wget.download(gff_url, out=gff_out)
            except Exception as e:
                print(f"Failed to download the file: {str(e)}")
        gff_dict = gff_parser( gff_out )
        os.remove(gff_out)

        # Downloading and parsing the eggNOG file of clusters size =1
        eggnog_url = rep_url+rep+'_eggNOG.tsv'
        eggnog_out = rep+'.eggNOG'
        if not os.path.exists(eggnog_out):
            try:
                wget.download(eggnog_url, out=eggnog_out)
            except Exception as g:
                print(f"Failed to download the file: {str(g)}")
        gff_dict = eggnog_parser( eggnog_out, gff_dict )
        os.remove(eggnog_out)

        # Generating output for cluster = 1
        if len(reps_clusters[rep]) == 1:
            output_writer( gff_dict, rep, 1, core_mgygs )

        # Parsing the accesory genes annotation
        else:
            acc_gff_dict = {}
            clstr_size = len(reps_clusters[rep])

            # Downloading and parsing the core genes list
            pan_url = url_prefix+'/species_catalogue/'+rep_prefix+'/'+rep+'/pan-genome/'
            core_tab_url = pan_url+'core_genes.txt'
            core_tab_out = rep+'_core.csv'
            if not os.path.exists(core_tab_out):
                try:
                    wget.download(core_tab_url, out=core_tab_out)
                except Exception as c:
                    print(f"Failed to download the file: {str(c)}")

            # Saving the core genes ids
            with open(core_tab_out, 'r') as input_file:
                for line in input_file:
                    gene_name = line.rstrip()
                    core_list.append(gene_name)
            os.remove(core_tab_out)

            # Downloading and parsing the pangenomic table
            acc_tab_url = pan_url+'gene_presence_absence.csv'
            acc_tab_out = rep+'_acc.csv'
            if not os.path.exists(acc_tab_out):
                try:
                    wget.download(acc_tab_url, out=acc_tab_out)
                except Exception as e:
                    print(f"Failed to download the file: {str(e)}")

            # Keep only one accessory gene per genome
            accesory_genes = {}
            relevant_members = []
            relevant_genes = []
            with open(acc_tab_out, 'r') as input_file:
                next(input_file)
                for line in input_file:
                    l_line = line.rstrip().split(",")
                    gene_key = l_line[0]
                    genomes_list = l_line[3:]
                    for member_gen in genomes_list:
                        if 'MGYG' in member_gen:
                            if ';' in member_gen:
                                for each_gene in member_gen.split(';'):
                                    if 'MGYG' in each_gene:
                                        prefix = member_gen.split('_')[0]
                            else:
                                prefix = member_gen.split('_')[0]
                                if prefix != rep:
                                    if not gene_key in accesory_genes:
                                        relevant_members.append(prefix)
                                        accesory_genes[gene_key]=member_gen
                                        relevant_genes.append(member_gen)
                                        if gene_key in core_list:
                                            core_mgygs.append(member_gen)
                                else:
                                    if gene_key in core_list:
                                        core_mgygs.append(member_gen)
            os.remove(acc_tab_out)

            # Downloading and saving coordinates and strand of selected accessory genes
            relevant_members = list(set(relevant_members))
            for member in relevant_members:
                all_url = url_prefix+'/all_genomes/'+rep_prefix+'/'+rep+'/genomes1/'
                mem_gff_url = all_url+member+'.gff.gz'
                mem_gff_out = member+'.gff.gz'
                if not os.path.exists(mem_gff_out):
                    try:
                        wget.download(mem_gff_url, out=mem_gff_out)
                    except Exception as e:
                        print(f"Failed to download the file: {str(e)}")

                acc_gff_dict = acc_gff_parser( mem_gff_out, relevant_genes, acc_gff_dict )
                os.remove(mem_gff_out)

            # Parsing the eggNOG annotation for accessory genes of relevant genomes
            eggnog_annot = './emapper_results/'+rep+'_out.emapper.annotations'
            if os.path.exists(eggnog_annot):
                acc_gff_dict = acc_eggnog_parser( eggnog_annot, accesory_genes, acc_gff_dict, pfam_desc )

            # Integrating representative and accessory annotations
            integrated_dict = {}
            for rep_gene, annot_data in gff_dict.items():
                integrated_dict[rep_gene] = annot_data
            for acc_gene, acc_annot_data in acc_gff_dict.items():
                if len(acc_annot_data) == 7:
                    integrated_dict[acc_gene] = acc_annot_data
            output_writer( integrated_dict, rep, clstr_size, core_mgygs )


def gff_parser( gff_file ):
    gff_dict = {}
    with open(gff_file, 'r') as input_file:
        for line in input_file:
            l_line = line.rstrip().split("\t")
            # Annotation lines have exactly 9 columns
            if len(l_line) == 9:
                (
                    contig,
                    seq_source,
                    seq_type,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    attr,
                ) = line.rstrip().split("\t")
                strand = strand.replace('+','1').replace('-','-1')
                if seq_type == 'CDS':
                    att_l = attr.split(";")
                    gene_id = att_l[0].replace('ID=','')
                    gff_dict[gene_id] = [contig, start, end, strand]
                    kegg_flag, pfam_flag = 0, 0
                    for attribute in att_l:
                        att_key,att_val = attribute.split('=')
                        if att_key == 'KEGG':
                            ko = att_val.replace('ko:','')
                            kegg_flag = 1
                        if att_key == 'Pfam':
                            pfam = att_val
                            pfam_flag = 1
                    if kegg_flag == 0:
                        ko = '-'
                    if pfam_flag == 0:
                        pfam = '-'
                    gff_dict[gene_id].append(ko)
                    gff_dict[gene_id].append(pfam)
    return( gff_dict )


def eggnog_parser( eggnog_out, gff_dict ):
    cazy_annot = {}
    with open(eggnog_out, 'r') as input_file:
        next( input_file )
        for line in input_file:
            l_line = line.rstrip().split("\t")
            gene_id = l_line[0]
            cazy = l_line[18]
            cazy_annot[gene_id] = cazy
    for gene in gff_dict:
        if gene in cazy_annot:
            gff_dict[gene].append(cazy_annot[gene])
        else:
            gff_dict[gene].append('-')
    return( gff_dict )


def acc_gff_parser( mem_gff_out, relevant_genes, acc_gff_dict ):
    with gzip.open(mem_gff_out, 'rt') as input_file:
        for line in input_file:
            l_line = line.rstrip().split("\t")
            # Annotation lines have exactly 9 columns
            if len(l_line) == 9:
                (
                    contig,
                    seq_source,
                    seq_type,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    attr,
                ) = line.rstrip().split("\t")
                strand = strand.replace('+','1').replace('-','-1')
                if seq_type == 'CDS':
                    att_l = attr.split(";")
                    gene_id = att_l[0].replace('ID=','')
                    if gene_id in relevant_genes:
                        acc_gff_dict[gene_id] = [contig, start, end, strand]
    return(acc_gff_dict)


def acc_eggnog_parser( eggnog_annot, accesory_genes, acc_gff_dict, pfam_desc ):
    rev_accesory_genes = {}
    for pan_gene, genome_gene in accesory_genes.items():
        rev_accesory_genes[genome_gene] = pan_gene
    kegg_annot = {}
    pfam_annot = {}
    cazy_annot = {}
    with open(eggnog_annot, 'r') as input_file:
        next( input_file )
        for line in input_file:
            l_line = line.rstrip().split("\t")
            gene_id = l_line[0]
            kegg_ko = l_line[11].replace('ko:','')
            kegg_annot[gene_id] = kegg_ko
            cazy = l_line[18]
            cazy_annot[gene_id] = cazy
            pfams_raw = l_line[20]
            if pfams_raw != '-':
                if ',' in pfams_raw:
                    pfam_ids_list = []
                    for each_pfam in pfams_raw.split(','):
                        if each_pfam in pfam_desc:
                            pfam_id = pfam_desc[each_pfam]
                            pfam_ids_list.append(pfam_id)
                    pfam_ids_list = list(set(pfam_ids_list))
                    if len(pfam_ids_list)>0:
                        pfams = ','.join(pfam_ids_list)
                    else:
                        pfams = '-'
                else:
                    if pfams_raw in pfam_desc:
                        pfams = pfam_desc[pfams_raw]
                    else:
                        pfams = '-'
            else:
                pfams = pfams_raw
            pfam_annot[gene_id] = pfams

    for gene in acc_gff_dict:
        pan_gene  = rev_accesory_genes[gene]
        if pan_gene in kegg_annot:
            acc_gff_dict[gene].append(kegg_annot[pan_gene])
        else:
            acc_gff_dict[gene].append('-')
        if pan_gene in pfam_annot:
            acc_gff_dict[gene].append(pfam_annot[pan_gene])
        else:
            acc_gff_dict[gene].append('-')
        if pan_gene in cazy_annot:
            acc_gff_dict[gene].append(cazy_annot[pan_gene])
        else:
            acc_gff_dict[gene].append('-')
    return(acc_gff_dict)



def output_writer( gff_dict, rep, clstr_size, core_mgygs ):
    with open(rep+'_clstr.tsv','w') as file_out:
        file_out.write('#cluster size = '+str(clstr_size)+'\n')
        file_out.write("\t".join([
            '#contig',
            'gene_id',
            'start',
            'end',
            'strand',
            'kegg',
            'pfam',
            'cazy',
            'core'
        ])+'\n')
        for gene in gff_dict:
            if clstr_size == 1:
                core = 'true'
            else:
                if gene in core_mgygs:
                    core = 'true'
                else:
                    core = 'false'
            file_out.write("\t".join([
                gff_dict[gene][0],
                gene,
                gff_dict[gene][1],
                gff_dict[gene][2],
                gff_dict[gene][3],
                gff_dict[gene][4],
                gff_dict[gene][5],
                gff_dict[gene][6],
                core
            ])+'\n')



def main():
    parser = argparse.ArgumentParser(
        description="This script find the accessory genes that needs eggNOG annotation. Provide the genomes catalogue metadata file and the url of the ftp site where the catalogue files are located"
    )
    parser.add_argument(
        "--metadata",
        type=str,
        help="Catalogue metadata file",
        required=True,
    )
    parser.add_argument(
        "--url_prefix",
        type=str,
        help="url of the ftp site where the catalogue files are located",
        required=True,
    )
    parser.add_argument(
        "--mode",
        type=str,
        help="Build profiles mode: pre or post. Option `pre` will generate the fasta files of the accesory genes to run eggNOG. Option `post` will process eggNOG outputs and download and parse the corresponding GFF files",
        required=True,
    )
    args = parser.parse_args()
    
    ### Calling functions
    ( pfam_desc ) = pfam_parser( '/hps/nobackup/rdf/metagenomics/service-team/ref-dbs/dram/Pfam-A.hmm.dat.gz' )
    ( reps_clusters ) = metadata_parser( args.metadata )

    if args.mode == 'pre':
        accessory_writer( reps_clusters, args.url_prefix )
    elif args.mode == 'post':
        annot_writer( reps_clusters, args.url_prefix, pfam_desc )
    else:
        print('Option '+args.mode+' is not valid. Provide a valid option on the mode argument: `pre` or `post`')
        exit()


if __name__ == "__main__":
    main()

