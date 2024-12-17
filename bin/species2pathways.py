#!/usr/bin/env python

import argparse


def relab_parser(relab_table):
    reps_list = []
    with open(relab_table) as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split("\t")
            rep_genome = l_line[0].split(";")[-1]
            reps_list.append(rep_genome)
    return reps_list


def pathways_finder(reps_list, kegg_comp_db, core_mode):
    species_values = {}
    all_pathways = []
    for rep_genome in reps_list:
        db_file = kegg_comp_db + "/" + rep_genome + "_clstr_kegg_comp.tsv"
        with open(db_file) as input_file:
            next(input_file)
            for line in input_file:
                module, pangenome, core = line.rstrip().split("\t")
                if core_mode == "pan":
                    completeness = pangenome
                elif core_mode == "core":
                    completeness = core
                else:
                    exit(
                        "The option "
                        + core_mode
                        + " is not allowed. Please use either `pan` or `core` strings"
                        + "\n"
                    )

                if float(completeness) > 0:
                    composite_key = (rep_genome, module)
                    species_values[composite_key] = completeness
                    all_pathways.append(module)

    all_pathways = list(set(all_pathways))

    return (all_pathways, species_values)


def species_writer(reps_list, all_pathways, species_values, output):
    genomes_names = [genome + "_clstr" for genome in reps_list]
    with open(output + "_species_kegg_modules_comp.tsv", "w") as output_file:
        header = ["module"] + genomes_names
        output_file.write("\t".join(header) + "\n")
        for pathway in all_pathways:
            to_print = []
            to_print.append(pathway)
            for species in reps_list:
                composite_key = (species, pathway)
                if composite_key in species_values:
                    to_print.append(species_values[composite_key])
                else:
                    to_print.append("0")
            output_file.write("\t".join(to_print) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="This script use the species prediction to generate functional tables from the pangenomic profiles"
    )
    parser.add_argument(
        "--kegg_comp_db",
        type=str,
        help="Path to the precomputed database of pangenomic kegg completeness profiles",
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
        help="Prefix to be used to name the output file",
        required=True,
    )
    args = parser.parse_args()

    ### Calling functions
    (reps_list) = relab_parser(args.relab)

    (all_pathways, species_values) = pathways_finder(
        reps_list, args.kegg_comp_db, args.core_mode
    )

    species_writer(reps_list, all_pathways, species_values, args.output)


if __name__ == "__main__":
    main()
