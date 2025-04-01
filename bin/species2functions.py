#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import argparse
import gzip


def pfam_parser(pfam_data):
    pfam_desc = {}
    with gzip.open(pfam_data, "rt", encoding="utf-8") as input_file:
        entry = []
        for line in input_file:
            line = line.strip()
            if line == "//" and entry:
                entry_dict = {}
                for entry_line in entry:
                    if entry_line.startswith("#=GF DE"):
                        desc_string = entry_line.split("#=GF DE   ")[1]
                        entry_dict["DE"] = desc_string
                    elif entry_line.startswith("#=GF AC"):
                        entry_dict["AC"] = entry_line.split(" ")[-1].split(".")[0]
                if "DE" in entry_dict and "AC" in entry_dict:
                    pfam_desc[entry_dict["AC"]] = entry_dict["DE"]
                entry = []
            else:
                entry.append(line)
    return pfam_desc


def dram_parser(dram_form):
    dram_desc = {}
    with open(dram_form) as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split("\t")
            gene_id = l_line[0]
            gene_description = l_line[1].replace('"', "")
            if gene_id in dram_desc:
                if gene_description not in dram_desc[gene_id]:
                    dram_desc[gene_id].append(gene_description)
            else:
                dram_desc[gene_id] = [gene_description]
    return dram_desc


def relab_parser(relab_table):
    taxonomy = {}
    reps_list = []
    with open(relab_table) as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split("\t")
            rep_genome = l_line[0].split(";")[-1]
            reps_list.append(rep_genome)
            lineage = l_line[0].split(";")
            lineage.pop(-1)
            lineage = [element for element in lineage if len(element) != 3]
            last_rank = lineage[-1]
            taxonomy[rep_genome] = last_rank
    return (reps_list, taxonomy)


def functions_finder_pan(reps_list, db_path):
    functions_dict = {}
    species_kos = {}
    species_pfams = {}
    per_gene_dict = {}
    gene_positions = {}
    all_kos = []
    all_pfams = []
    for rep_genome in reps_list:
        db_file = db_path + "/" + rep_genome + "_clstr.tsv"
        positions = {}
        species_kos[rep_genome] = []
        species_pfams[rep_genome] = []
        per_gene_dict[rep_genome] = []
        pan_kos = []
        pan_pfams = []
        with open(db_file) as input_file:
            next(input_file)
            next(input_file)
            for line in input_file:
                per_gene_dict[rep_genome].append(line.rstrip())
                (
                    contig,
                    gene_id,
                    start,
                    end,
                    strand,
                    kegg,
                    cazy,
                    pfam,
                    core,
                ) = line.rstrip().split("\t")

                if contig not in positions:
                    positions[contig] = {}
                    pos = 1
                else:
                    pos += 1
                gene_positions[gene_id] = str(pos)

                if kegg != "-":
                    if "," in kegg:
                        kegg_list = kegg.split(",")
                    else:
                        kegg_list = [kegg]
                    pan_kos = pan_kos + kegg_list
                    for current_ko in kegg_list:
                        if current_ko not in species_kos[rep_genome]:
                            species_kos[rep_genome].append(current_ko)

                if pfam != "-":
                    if "," in pfam:
                        pfam_list = pfam.split(",")
                    else:
                        pfam_list = [pfam]
                    pan_pfams = pan_pfams + pfam_list
                    for current_pfam in pfam_list:
                        if current_pfam not in species_pfams[rep_genome]:
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

    return (
        functions_dict,
        all_kos,
        all_pfams,
        species_kos,
        species_pfams,
        per_gene_dict,
        gene_positions,
    )


def functions_finder_core(reps_list, db_path):
    functions_dict = {}
    species_kos = {}
    species_pfams = {}
    per_gene_dict = {}
    gene_positions = {}
    all_kos = []
    all_pfams = []
    for rep_genome in reps_list:
        db_file = db_path + "/" + rep_genome + "_clstr.tsv"
        positions = {}
        species_kos[rep_genome] = []
        species_pfams[rep_genome] = []
        per_gene_dict[rep_genome] = []
        pan_kos = []
        pan_pfams = []
        with open(db_file) as input_file:
            next(input_file)
            next(input_file)
            for line in input_file:
                (
                    contig,
                    gene_id,
                    start,
                    end,
                    strand,
                    kegg,
                    cazy,
                    pfam,
                    core,
                ) = line.rstrip().split("\t")

                if contig not in positions:
                    positions[contig] = {}
                    pos = 1
                else:
                    pos += 1
                gene_positions[gene_id] = str(pos)

                if core == "true":
                    per_gene_dict[rep_genome].append(line.rstrip())
                    if kegg != "-":
                        if "," in kegg:
                            kegg_list = kegg.split(",")
                        else:
                            kegg_list = [kegg]
                        pan_kos = pan_kos + kegg_list
                        for current_ko in kegg_list:
                            if current_ko not in species_kos[rep_genome]:
                                species_kos[rep_genome].append(current_ko)

                    if pfam != "-":
                        if "," in pfam:
                            pfam_list = pfam.split(",")
                        else:
                            pfam_list = [pfam]
                        pan_pfams = pan_pfams + pfam_list
                        for current_pfam in pfam_list:
                            if current_pfam not in species_pfams[rep_genome]:
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

    return (
        functions_dict,
        all_kos,
        all_pfams,
        species_kos,
        species_pfams,
        per_gene_dict,
        gene_positions,
    )


def community_writer(functions_dict, all_kos, all_pfams, output):
    with open(output + "_community_kegg.tsv", "w") as output_file:
        output_file.write("ko_id\t" + output + "\n")
        for ko in all_kos:
            output_file.write(ko + "\t" + str(functions_dict[ko]) + "\n")

    with open(output + "_community_pfams.tsv", "w") as output_file:
        output_file.write("pfam_id\t" + output + "\n")
        for pfam in all_pfams:
            output_file.write(pfam + "\t" + str(functions_dict[pfam]) + "\n")


def species_writer(reps_list, all_kos, all_pfams, species_kos, species_pfams, output):
    genomes_names = [genome + "_clstr" for genome in reps_list]

    with open(output + "_species_kegg.tsv", "w") as output_file:
        header = ["ko_id"] + genomes_names
        output_file.write("\t".join(header) + "\n")
        for ko in all_kos:
            to_print = []
            to_print.append(ko)
            for species in reps_list:
                if ko in species_kos[species]:
                    to_print.append("1")
                else:
                    to_print.append("0")
            output_file.write("\t".join(to_print) + "\n")

    with open(output + "_species_pfams.tsv", "w") as output_file:
        header = ["pfam_id"] + genomes_names
        output_file.write("\t".join(header) + "\n")
        for pfam in all_pfams:
            to_print = []
            to_print.append(pfam)
            for species in reps_list:
                if pfam in species_pfams[species]:
                    to_print.append("1")
                else:
                    to_print.append("0")
            output_file.write("\t".join(to_print) + "\n")


def dram_writer(per_gene_dict, gene_positions, taxonomy, pfam_desc, dram_desc, output):
    species_annot = {}
    dram_header = [
        "",
        "fasta",
        "scaffold",
        "gene_position",
        "kegg_id",
        "kegg_hit",
        "pfam_hits",
        "cazy_hits",
        "bin_taxonomy",
    ]
    counter = 0

    with open(output + "_species_dram.tsv", "w") as output_sp, open(
        output + "_community_dram.tsv", "w"
    ) as output_comm:
        output_sp.write("\t".join(dram_header) + "\n")

        # Parsing the genes dictionary. We use the species_clstr id instead of fasta
        for species_clstr in per_gene_dict:
            cazy_hits, pfam_hits, kegg_ids, kegg_hits = [], [], [], []
            counter += 1

            species_annot[species_clstr] = {}
            species_annot[species_clstr]["fasta"] = (
                    species_clstr + "_clstr: " + taxonomy[species_clstr]
            )
            species_annot[species_clstr]["scaffold"] = species_clstr + '_' + str(counter)
            species_annot[species_clstr]["gene_position"] = species_clstr + '_clstr'
            species_annot[species_clstr]["bin_taxonomy"] = taxonomy[species_clstr]

            for gene_line in per_gene_dict[species_clstr]:
                # Populating handy info
                (
                    contig,
                    gene_id,
                    start,
                    end,
                    strand,
                    kegg,
                    cazy,
                    pfam,
                    core,
                ) = gene_line.split("\t")

                # Processing cazy, pfam, and kegg descriptions. If no description,
                # then the annotation is discarded to avoid passing depricated annotation to DRAM
                #rank = "E"
                if cazy != '-':
                    cazy_desc_list = []
                    for cazy_acc in cazy.split(","):
                        if cazy_acc in dram_desc:
                            for cazy_desc in dram_desc[cazy_acc]:
                                if any(
                                    [cazy_desc.endswith(";"), cazy_desc.endswith(".")]
                                ):
                                    cazy_desc = cazy_desc[:-1]
                                cazy_desc_list.append(cazy_desc)
                            last_element = cazy_desc_list.pop(-1)
                            last_element = last_element + " [" + cazy_acc + "]"
                            if last_element not in cazy_hits:
                                cazy_hits.append(last_element)

                if pfam != '-':
                    for pfam_id in pfam.split(","):
                        if pfam_id in pfam_desc:
                            pfam_full_desc = pfam_desc[pfam_id] + " [" + pfam_id + "]"
                            if pfam_full_desc not in pfam_hits:
                                pfam_hits.append(pfam_full_desc)

                if kegg != "-":
                    for ko in kegg.split(","):
                        if ko in dram_desc:
                            if ko not in kegg_ids:
                                kegg_ids.append(ko)
                            for ko_desc in dram_desc[ko]:
                                if ko_desc not in kegg_hits:
                                    kegg_hits.append(ko_desc)

                # Aggregating annotation at genome level for species output
                species_annot[species_clstr].setdefault("cazy_hits", []).extend(cazy_hits)
                species_annot[species_clstr].setdefault("pfam_hits", []).extend(pfam_hits)
                species_annot[species_clstr].setdefault("kegg_id", []).extend(kegg_ids)
                species_annot[species_clstr].setdefault("kegg_hit", []).extend(kegg_hits)


        # Writing output at assembly level and aggreagting at sample level sample
        sample_annot = {}
        for genome in species_annot:
            to_print = []
            to_print.append(genome)
            for header_key in dram_header[1:]:
                value = species_annot[genome][header_key]
                if isinstance(value, list):
                    value = set(value)
                    sample_annot.setdefault(header_key, []).extend(value)
                    value = "; ".join(value)
                to_print.append(value)
            output_sp.write("\t".join(to_print) + "\n")

        # Aggregating annotation at sample level for community output
        sample_annot["fasta"] = "community: " + output
        sample_annot["bin_taxonomy"] = "community_" + output
        sample_annot["scaffold"] = "community_0"
        sample_annot["gene_position"] = 'microbial_community'

        to_print = []
        to_print.append("assembly")
        for header_key in dram_header[1:]:
            value = sample_annot[header_key]
            if isinstance(value, list):
                value = "; ".join(set(value))
            to_print.append(value)
        output_comm.write("\t".join(to_print) + "\n")



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
        "--dram_out",
        action="store_true",
        help="Generate dram files for DRAM_distill. Default = false",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Prefix to be used to name the output files. This string will also be used in headers",
        required=True,
    )
    args = parser.parse_args()

    pfam_db = args.external_db + "/Pfam-A.hmm.dat.gz"

    ### Calling functions
    (pfam_desc) = pfam_parser(pfam_db)
    (reps_list, taxonomy) = relab_parser(args.relab)

    if args.core_mode == "pan":
        (
            functions_dict,
            all_kos,
            all_pfams,
            species_kos,
            species_pfams,
            per_gene_dict,
            gene_positions,
        ) = functions_finder_pan(reps_list, args.pangenome_db)
    elif args.core_mode == "core":
        (
            functions_dict,
            all_kos,
            all_pfams,
            species_kos,
            species_pfams,
            per_gene_dict,
            gene_positions,
        ) = functions_finder_core(reps_list, args.pangenome_db)
    else:
        exit(
            "The option "
            + args.core_mode
            + " is not allowed. Chose eaither `pan` or `core`"
        )

    community_writer(functions_dict, all_kos, all_pfams, args.output)
    species_writer(
        reps_list, all_kos, all_pfams, species_kos, species_pfams, args.output
    )

    if args.dram_out:
        dram_form = args.external_db + "/genome_summary_form.tsv"
        (dram_desc) = dram_parser(dram_form)
        dram_writer(
            per_gene_dict, gene_positions, taxonomy, pfam_desc, dram_desc, args.output
        )


if __name__ == "__main__":
    main()
