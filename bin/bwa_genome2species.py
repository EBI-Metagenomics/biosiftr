#!/usr/bin/env python

import argparse
import os.path
import sys
from Bio import SeqIO

##### This script transforms BWA genomes relative abundance into species relative abundance
##### Alejandra Escobar, EMBL-EBI
##### Dec 18, 2023


def metadata_parser(catalogue_metadata):
    ref_spec_genome = {}
    with open(catalogue_metadata, "r") as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split("\t")
            genome = l_line[13]
            lineage = l_line[14] + ";" + genome
            if not genome in ref_spec_genome:
                ref_spec_genome[genome] = lineage

    return ref_spec_genome


def aggregate_species(genomes_relab, ref_spec_genome, out_name):
    species_reads = {}
    total_reads = 0
    with open(genomes_relab, "r") as input_file:
        next(input_file)
        for line in input_file:
            (
                reference_id,
                reference_len,
                reads_count,
                base_cov,
                rel_abun,
            ) = line.rstrip().split("\t")

            total_reads = total_reads + int(reads_count)

            if reference_id in ref_spec_genome:
                lineage = ref_spec_genome[reference_id]
                if lineage in species_reads:
                    species_reads[lineage] = species_reads[lineage] + int(reads_count)
                else:
                    species_reads[lineage] = int(reads_count)

            else:
                print("No species found in metadata file for " + reference_id)

    with open(out_name, "w") as file_out:
        prefix_name = out_name.replace("_species.tsv", "")
        file_out.write("\t".join(["lineage", prefix_name]) + "\n")
        for species in species_reads:
            relab = float(species_reads[species]) / float(total_reads)
            file_out.write("\t".join([species.replace(" ", "_"), str(relab)]) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="This script transforms BWA genomes relative abundance into species relative abundance. Provide the genomes relative abundance file of unique mapping reads and the genomes catalogue metadata file"
    )
    parser.add_argument(
        "--genomes_relab",
        type=str,
        help="Relative abundance table generated from unique HQ mapping reads and cov thres >= 0.1",
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
        help="Name of the output file. Default = species_relab.tsv",
        required=False,
    )
    args = parser.parse_args()

    if args.output:
        out_name = args.output
    else:
        out_name = "species_relab.tsv"

    ### Calling functions
    (ref_spec_genome) = metadata_parser(args.metadata)
    aggregate_species(args.genomes_relab, ref_spec_genome, out_name)


if __name__ == "__main__":
    main()
