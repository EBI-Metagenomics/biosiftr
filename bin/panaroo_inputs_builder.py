#!/usr/bin/env python

import argparse
import os.path
import sys
import wget
import gzip
import shutil
import subprocess
from Bio import SeqIO


##### This script prepare the inputs to launch panaroo on the human-gut catalogue v2.0
##### Alejandra Escobar, EMBL-EBI
##### June 21, 2024


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


def unzip_gz_files(src_dir, dest_dir):
    for root, _, files in os.walk(src_dir):
        for file_name in files:
            if file_name.endswith(".gz"):
                full_file_name = os.path.join(root, file_name)
                relative_path = os.path.relpath(full_file_name, src_dir)
                dest_path = os.path.join(dest_dir, relative_path[:-3])

                # Ensure the destination directory exists
                os.makedirs(os.path.dirname(dest_path), exist_ok=True)

                with gzip.open(full_file_name, "rb") as f_in:
                    with open(dest_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)


def gff_gather(reps_clusters, loc_prefix):
    cluster_counter = 0
    total_clusters = len(reps_clusters)
    for rep in reps_clusters:
        cluster_counter += 1
        print(
            "Processing cluster number "
            + str(cluster_counter)
            + " out of "
            + str(total_clusters)
        )

        if len(reps_clusters[rep]) > 1:
            dest_dir = "rep_" + rep

            if rep.endswith(".1"):
                rep_prefix = rep[:-4]
            else:
                rep_prefix = rep[:-2]

            # Gathering the gff files
            src_dir = (
                loc_prefix + "/all_genomes/" + rep_prefix + "/" + rep + "/genomes1/"
            )

            print("Unzipping " + str(len(reps_clusters[rep])) + " files for " + rep)
            unzip_gz_files(src_dir, dest_dir)


def main():
    parser = argparse.ArgumentParser(
        description="This script prepare the inputs to launch panaroo on the human-gut catalogue v2.0"
    )
    parser.add_argument(
        "--metadata",
        type=str,
        help="Catalogue metadata file location",
        required=True,
    )
    args = parser.parse_args()

    ### Calling functions
    loc_prefix = args.metadata.replace("genomes-all_metadata.tsv", "")
    (reps_clusters) = metadata_parser(args.metadata)
    gff_gather(reps_clusters, loc_prefix)


if __name__ == "__main__":
    main()
