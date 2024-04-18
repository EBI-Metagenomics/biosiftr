#!/usr/bin/env python

import argparse
import os.path
import sys
import wget
import gzip
from Bio import SeqIO

##### This script find the accessory genes that needs eggNOG annotation in codon catalogues
##### Alejandra Escobar, EMBL-EBI
##### March 27, 2024


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


def metadata_parser(catalogue_metadata):
    reps_clusters = {}
    with open(catalogue_metadata, "r") as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split("\t")
            clstr_member = l_line[0]
            rep_genome = l_line[13]
            if rep_genome in reps_clusters:
                reps_clusters[rep_genome].append(clstr_member)
            else:
                reps_clusters[rep_genome] = [clstr_member]
    return reps_clusters


def accessory_writer(reps_clusters, loc_prefix):
    for rep in reps_clusters:
        if len(reps_clusters[rep]) > 1:
            if rep.endswith('.1'):
                rep_prefix = rep[:-4]
            else:
                rep_prefix = rep[:-2]
            pan_loc = (
                loc_prefix
                + "/species_catalogue/"
                + rep_prefix
                + "/"
                + rep
                + "/pan-genome/"
            )

            # Parsing the presence/absence tab
            r_tab_loc = pan_loc + "gene_presence_absence.Rtab"
            accessory_genes = []
            with open(r_tab_loc, "r") as input_file:
                header = input_file.readline().strip().split("\t")
                index = header.index(rep)
                for line in input_file:
                    l_line = line.rstrip().split("\t")
                    if l_line[index] == "0":
                        accessory_genes.append(l_line[0])

            # Parsing the fasta file of genes
            genes_loc = pan_loc + "pan-genome.fna"
            with open(rep + "_accessory.fasta", "w") as fasta_out:
                for record in SeqIO.parse(genes_loc, "fasta"):
                    seq_id = str(record.id)
                    if seq_id in accessory_genes:
                        fasta_out.write(">" + seq_id + "\n")
                        fasta_out.write(str(record.seq).upper() + "\n")


def annot_writer(reps_clusters, loc_prefix, pfam_desc):
    for rep in reps_clusters:
        core_list, core_mgygs = [], []
        if rep.endswith('.1'):
            rep_prefix = rep[:-4]
        else:
            rep_prefix = rep[:-2]
        rep_loc = (
            loc_prefix + "/species_catalogue/" + rep_prefix + "/" + rep + "/genome/"
        )

        # Parsing the gff file of clusters size = 1
        gff_loc = rep_loc + rep + ".gff"
        gff_dict = gff_parser(gff_loc)

        # Parsing the eggNOG file of clusters size = 1
        eggnog_loc = rep_loc + rep + "_eggNOG.tsv"
        gff_dict = eggnog_parser(eggnog_loc, gff_dict, pfam_desc)

        # Generating output for cluster = 1
        if len(reps_clusters[rep]) == 1:
            output_writer(gff_dict, rep, 1, core_mgygs)

        # Parsing the accesory genes annotation
        else:
            acc_gff_dict = {}
            clstr_size = len(reps_clusters[rep])

            # Parsing the core genes list
            pan_loc = (
                loc_prefix
                + "/species_catalogue/"
                + rep_prefix
                + "/"
                + rep
                + "/pan-genome/"
            )
            core_tab_loc = pan_loc + "core_genes.txt"

            # Saving the core genes ids
            with open(core_tab_loc, "r") as input_file:
                for line in input_file:
                    gene_name = line.rstrip()
                    core_list.append(gene_name)

            # Parsing the pangenomic table
            acc_tab_loc = pan_loc + "gene_presence_absence.csv"

            # Keep only one accessory gene per genome
            accesory_genes = {}
            relevant_members = []
            relevant_genes = []
            with open(acc_tab_loc, "r") as input_file:
                next(input_file)
                for line in input_file:
                    l_line = line.rstrip().split(",")
                    gene_key = l_line[0]
                    genomes_list = l_line[3:]
                    for member_gen in genomes_list:
                        if "MGYG" in member_gen:
                            if ";" in member_gen:
                                for each_gene in member_gen.split(";"):
                                    if "MGYG" in each_gene:
                                        prefix = member_gen.split("_")[0]
                            else:
                                prefix = member_gen.split("_")[0]
                                if prefix != rep:
                                    if not gene_key in accesory_genes:
                                        relevant_members.append(prefix)
                                        accesory_genes[gene_key] = member_gen
                                        relevant_genes.append(member_gen)
                                        if gene_key in core_list:
                                            core_mgygs.append(member_gen)
                                else:
                                    if gene_key in core_list:
                                        core_mgygs.append(member_gen)

            # Saving coordinates and strand of selected accessory genes
            relevant_members = list(set(relevant_members))
            for member in relevant_members:
                all_loc = (
                    loc_prefix + "/all_genomes/" + rep_prefix + "/" + rep + "/genomes1/"
                )
                mem_gff_loc = all_loc + member + ".gff.gz"
                acc_gff_dict = acc_gff_parser(mem_gff_loc, relevant_genes, acc_gff_dict)

            # Parsing the eggNOG annotation for accessory genes of relevant genomes
            eggnog_annot = "./emapper_results/" + rep + "_out.emapper.annotations"
            if os.path.exists(eggnog_annot):
                acc_gff_dict = acc_eggnog_parser(
                    eggnog_annot, accesory_genes, acc_gff_dict, pfam_desc
                )

            # Integrating representative and accessory annotations
            integrated_dict = {}
            for rep_gene, annot_data in gff_dict.items():
                integrated_dict[rep_gene] = annot_data
            for acc_gene, acc_annot_data in acc_gff_dict.items():
                if len(acc_annot_data) == 7:
                    integrated_dict[acc_gene] = acc_annot_data
            output_writer(integrated_dict, rep, clstr_size, core_mgygs)


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
    with gzip.open(mem_gff_out, "rt") as input_file:
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


def acc_eggnog_parser(eggnog_annot, accesory_genes, acc_gff_dict, pfam_desc):
    rev_accesory_genes = {}
    for pan_gene, genome_gene in accesory_genes.items():
        rev_accesory_genes[genome_gene] = pan_gene
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
        pan_gene = rev_accesory_genes[gene]
        if pan_gene in kegg_annot:
            acc_gff_dict[gene].append(kegg_annot[pan_gene])
        else:
            acc_gff_dict[gene].append("-")
        if pan_gene in cazy_annot:
            acc_gff_dict[gene].append(cazy_annot[pan_gene])
        else:
            acc_gff_dict[gene].append("-")
        if pan_gene in pfam_annot:
            acc_gff_dict[gene].append(pfam_annot[pan_gene])
        else:
            acc_gff_dict[gene].append("-")
    return acc_gff_dict


def output_writer(gff_dict, rep, clstr_size, core_mgygs):
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
                if gene in core_mgygs:
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
        description="This script find the accessory genes that needs eggNOG annotation. Provide the genomes catalogue metadata file and the url of the ftp site where the catalogue files are located"
    )
    parser.add_argument(
        "--metadata",
        type=str,
        help="Catalogue metadata file location",
        required=True,
    )
    parser.add_argument(
        "--mode",
        type=str,
        help="Build profiles mode: pre or post. Option `pre` will generate the fasta files of the accesory genes to run eggNOG. Option `post` will process eggNOG outputs and download and parse the corresponding GFF files",
        required=True,
    )
    parser.add_argument(
        "--pfam_dat",
        type=str,
        help="The Pfam-A.hmm.dat.gz file location",
        required=True,
    )
    args = parser.parse_args()

    ### Calling functions
    (pfam_desc) = pfam_parser(args.pfam_dat)
    (reps_clusters) = metadata_parser(args.metadata)

    loc_prefix = args.metadata.replace("genomes-all_metadata.tsv", "")

    if args.mode == "pre":
        accessory_writer(reps_clusters, loc_prefix)
    elif args.mode == "post":
        annot_writer(reps_clusters, loc_prefix, pfam_desc)
    else:
        exit(
            "Option "
            + args.mode
            + " is not valid. Provide a valid option on the mode argument: `pre` or `post`"
        )


if __name__ == "__main__":
    main()
