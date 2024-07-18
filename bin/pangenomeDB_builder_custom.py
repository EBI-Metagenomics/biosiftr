#!/usr/bin/env python

import argparse
import os.path
import sys
import wget
import gzip
from Bio import SeqIO

##### This script find the accessory genes that needs eggNOG annotation for custom genome catalogues
##### Alejandra Escobar, EMBL-EBI
##### July 10, 2024


def pfam_parser(pfam_data):
    pfam_desc = {}
    with gzip.open(pfam_data, "rt", encoding="utf-8") as input_file:
        entry = []
        for line in input_file:
            line = line.strip()
            if line == "//" and entry:
                entry_dict = {}
                for entry_line in entry:
                    if entry_line.startswith("#=GF ID"):
                        entry_dict["ID"] = entry_line.split()[-1]
                    elif entry_line.startswith("#=GF AC"):
                        entry_dict["AC"] = entry_line.split()[-1].split(".")[0]
                if "ID" in entry_dict and "AC" in entry_dict:
                    pfam_desc[entry_dict["ID"]] = entry_dict["AC"]
                entry = []
            else:
                entry.append(line)
    return pfam_desc


def metadata_parser( drep_clstrs, derep_genomes ):
    clusters = {}
    with open(drep_clstrs, "r") as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split(",")
            clstr_member = l_line[0].replace('.fa','')
            clstr_num = 'clstr_'+l_line[5]
            if clstr_num in clusters:
                clusters[clstr_num].append(clstr_member)
            else:
                clusters[clstr_num] = [clstr_member]

    reps_clusters = {}
    rep_assemblies_list = os.listdir(derep_genomes)
    rep_genomes_list = [f.replace('.fa', '') for f in rep_assemblies_list]
    for clstr in clusters:
        members = clusters[clstr]
        for genome in members:
            if genome in rep_genomes_list:
                reps_clusters[genome] = members

    return reps_clusters


def annot_writer(reps_clusters, prokka_path, panaroo_path, pfam_desc):
    for rep in reps_clusters:
        core_list = []
        rep_loc = prokka_path + '/' + rep + '_prokka/'

        # Parsing the gff file of clusters size = 1
        gff_loc = rep_loc + rep + ".gff"
        gff_dict = gff_parser(gff_loc)

        # Parsing the eggNOG file of clusters size = 1
        eggnog_loc = rep_loc + rep + "_eggnog.emapper.annotations"
        gff_dict = eggnog_parser(eggnog_loc, gff_dict, pfam_desc)

        # Generating output for clusters size = 1
        if len(reps_clusters[rep]) == 1:
            output_writer(gff_dict, rep, 1, core_list)

        # Parsing pangenome of clusters size > 1
        else:
            relevant_members = [rep]
            clstr_size = len(reps_clusters[rep])

            # Parsing the core genes list
            pan_loc = (
                panaroo_path
                + "/rep_"
                + rep
                + "_panaroo/"
            )
            core_tab_loc = pan_loc + "core_genes.txt"

            # Saving the core genes ids
            with open(core_tab_loc, "r") as input_file:
                for line in input_file:
                    gene_name = line.rstrip()
                    core_list.append(gene_name)

            # Parsing the pangenomic table to keep only one accessory gene per genome
            acc_tab_loc = pan_loc + "gene_presence_absence.csv"
            pan_genes = {}
            with open(acc_tab_loc, "r") as input_file:
                next(input_file)
                for line in input_file:
                    l_line = line.rstrip().split(",")
                    gene_key = l_line[0]
                    genomes_list = l_line[3:]
                    pan_genes[gene_key] = []
                    for member_gen in genomes_list:
                        if len(member_gen) > 0:
                            if ";" in member_gen:
                                paralogs = [p for p in member_gen.split(";") if 'refound' not in p]
                                if len(paralogs) > 0:
                                    first_copy = paralogs[0]
                                    pan_genes[gene_key].append(first_copy)
                            elif not "refound" in member_gen:
                               pan_genes[gene_key].append(member_gen) 

            # Giving priority to genes in the representative genome
            relevant_genes = {}
            core_names = []
            for gene_key in pan_genes:
                rep_gene = [gene for gene in pan_genes[gene_key] if rep in gene]
                if len(rep_gene) == 0:
                    selected_gene = pan_genes[gene_key][0]
                    prefix = selected_gene.split('_')
                    prefix.pop(-1)
                    prefix = '_'.join(prefix)
                    relevant_members.append(prefix)
                else:
                    selected_gene = rep_gene[0]
                relevant_genes[selected_gene] = gene_key
                if gene_key in core_list:
                    core_names.append(selected_gene)

            # Saving coordinates and strand of selected genes
            relevant_members = list(set(relevant_members))
            
            acc_gff_dict = {}
            for member in relevant_members:
                annot_loc = prokka_path + '/' + member + '_prokka/'

                # Parsing the gff file of relevant members including the rep genome
                mem_gff_loc = annot_loc + member + ".gff"
                acc_gff_dict = acc_gff_parser(mem_gff_loc, relevant_genes, acc_gff_dict)

                # Parsing the eggNOG annotation of relevant genes
                eggnog_annot = annot_loc + member + "_eggnog.emapper.annotations"
                if os.path.exists(eggnog_annot):
                    acc_gff_dict = acc_eggnog_parser(eggnog_annot, relevant_genes, acc_gff_dict, pfam_desc)

                

            
            '''
            # Integrating representative and accessory annotations
            integrated_dict = {}
            for rep_gene, annot_data in gff_dict.items():
                integrated_dict[rep_gene] = annot_data
            for acc_gene, acc_annot_data in acc_gff_dict.items():
                if len(acc_annot_data) == 7:
                    integrated_dict[acc_gene] = acc_annot_data

            '''

            output_writer(acc_gff_dict, rep, clstr_size, core_names)
            
            


def gff_parser(gff_file):
    gff_dict = {}
    with open(gff_file, "r") as input_file:
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
                strand = strand.replace("+", "1").replace("-", "-1")
                if seq_type == "CDS":
                    att_l = attr.split(";")
                    gene_id = att_l[0].replace("ID=", "")
                    gff_dict[gene_id] = [contig, start, end, strand]
    return gff_dict


def eggnog_parser(eggnog_out, gff_dict, pfam_desc):
    ko_annot, cazy_annot, pfam_annot = {}, {}, {}
    with open(eggnog_out, "r") as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split("\t")
            gene_id = l_line[0]
            evalue = float(l_line[2])
            if evalue <= 1e-10:
                kegg_ko = l_line[11].replace("ko:", "")
                ko_annot[gene_id] = kegg_ko
                cazys = l_line[18]
                cazy_annot[gene_id] = cazys
                pfams_raw = l_line[20]
                if pfams_raw != "-":
                    if "," in pfams_raw:
                        pfam_ids_list = []
                        for each_pfam in pfams_raw.split(","):
                            if each_pfam in pfam_desc:
                                pfam_id = pfam_desc[each_pfam]
                                pfam_ids_list.append(pfam_id)
                        pfam_ids_list = list(set(pfam_ids_list))
                        if len(pfam_ids_list) > 0:
                            pfams = ",".join(pfam_ids_list)
                        else:
                            pfams = "-"
                    else:
                        if pfams_raw in pfam_desc:
                            pfams = pfam_desc[pfams_raw]
                        else:
                            pfams = "-"
                else:
                    pfams = pfams_raw
                pfam_annot[gene_id] = pfams

    for gene in gff_dict:
        if gene in ko_annot:
            gff_dict[gene].append(ko_annot[gene])
        else:
            gff_dict[gene].append("-")
        if gene in cazy_annot:
            gff_dict[gene].append(cazy_annot[gene])
        else:
            gff_dict[gene].append("-")
        if gene in pfam_annot:
            gff_dict[gene].append(pfam_annot[gene])
        else:
            gff_dict[gene].append("-")

    return gff_dict


def acc_gff_parser(mem_gff_out, relevant_genes, acc_gff_dict):
    with open(mem_gff_out, "r") as input_file:
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
                strand = strand.replace("+", "1").replace("-", "-1")
                if seq_type == "CDS":
                    att_l = attr.split(";")
                    gene_id = att_l[0].replace("ID=", "")
                    if gene_id in relevant_genes:
                        acc_gff_dict[gene_id] = [contig, start, end, strand]
    return acc_gff_dict


def acc_eggnog_parser(eggnog_annot, relevant_genes, acc_gff_dict, pfam_desc):
    kegg_annot, pfam_annot, cazy_annot = {}, {}, {}

    with open(eggnog_annot, "r") as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split("\t")
            gene_id = l_line[0]
            evalue = float(l_line[2])
            if evalue <= 1e-10:
                kegg_ko = l_line[11].replace("ko:", "")
                kegg_annot[gene_id] = kegg_ko
                cazy = l_line[18]
                cazy_annot[gene_id] = cazy
                pfams_raw = l_line[20]
                if pfams_raw != "-":
                    if "," in pfams_raw:
                        pfam_ids_list = []
                        for each_pfam in pfams_raw.split(","):
                            if each_pfam in pfam_desc:
                                pfam_id = pfam_desc[each_pfam]
                                pfam_ids_list.append(pfam_id)
                        pfam_ids_list = list(set(pfam_ids_list))
                        if len(pfam_ids_list) > 0:
                            pfams = ",".join(pfam_ids_list)
                        else:
                            pfams = "-"
                    else:
                        if pfams_raw in pfam_desc:
                            pfams = pfam_desc[pfams_raw]
                        else:
                            pfams = "-"
                else:
                    pfams = pfams_raw
                pfam_annot[gene_id] = pfams

    for gene in acc_gff_dict:
        if gene in relevant_genes:
            if gene in kegg_annot:
                acc_gff_dict[gene].append(kegg_annot[gene])
            else:
                acc_gff_dict[gene].append("-")
            if gene in cazy_annot:
                acc_gff_dict[gene].append(cazy_annot[gene])
            else:
                acc_gff_dict[gene].append("-")
            if gene in pfam_annot:
                acc_gff_dict[gene].append(pfam_annot[gene])
            else:
                acc_gff_dict[gene].append("-")

    return acc_gff_dict



def output_writer(gff_dict, rep, clstr_size, core_names):
    with open(rep + "_clstr.tsv", "w") as file_out:
        file_out.write("#cluster size = " + str(clstr_size) + "\n")
        file_out.write(
            "\t".join(
                [
                    "#contig",
                    "gene_id",
                    "start",
                    "end",
                    "strand",
                    "kegg",
                    "cazy",
                    "pfam",
                    "core",
                ]
            )
            + "\n"
        )
        for gene in gff_dict:
            if clstr_size == 1:
                core = "true"
            else:
                if gene in core_names:
                    core = "true"
                else:
                    core = "false"
            file_out.write(
                "\t".join(
                    [
                        gff_dict[gene][0],
                        gene,
                        gff_dict[gene][1],
                        gff_dict[gene][2],
                        gff_dict[gene][3],
                        gff_dict[gene][4],
                        gff_dict[gene][5],
                        gff_dict[gene][6],
                        core,
                    ]
                )
                + "\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="This script find the accessory genes that needs eggNOG annotation for custom genome catalogues"
    )
    parser.add_argument(
        "--drep_clstrs",
        type=str,
        help="Data table generated by drep: drep_output/data_tables/Cdb.csv",
        required=True,
    )
    parser.add_argument(
        "--derep_genomes",
        type=str,
        help="Path to the dereplicated genomes: drep_output/dereplicated_genomes/",
        required=True,
    )
    parser.add_argument(
        "--pfam_dat",
        type=str,
        help="The Pfam-A.hmm.dat.gz file location",
        required=True,
    )
    parser.add_argument(
        "--panaroo_path",
        type=str,
        help="The path to panaroo results",
        required=True,
    )
    parser.add_argument(
        "--prokka_path",
        type=str,
        help="The path to prokka results. EggNOG annotations should exists there",
        required=True,
    )
    args = parser.parse_args()

    ### Calling functions
    pfam_desc = pfam_parser( args.pfam_dat )
    reps_clusters = metadata_parser( args.drep_clstrs, args.derep_genomes )
    annot_writer(reps_clusters, args.prokka_path, args.panaroo_path, pfam_desc)


if __name__ == "__main__":
    main()
