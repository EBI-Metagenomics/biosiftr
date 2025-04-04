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


def matrix_parser(matrix_file, features_dict, all_features, all_samples):
    with open(matrix_file) as input_file:
        samples_list = input_file.readline().strip().split("\t")
        samples_list.pop(0)
        all_samples = all_samples + samples_list
        for line in input_file:
            l_line = line.rstrip().split("\t")
            feature = l_line.pop(0)
            all_features.append(feature)
            index = 0
            for element in l_line:
                sample_name = samples_list[index]
                composite_key = (sample_name, feature)
                index += 1
                if composite_key in features_dict:
                    print(
                        "Sample "
                        + sample_name
                        + " was found multiple times, it cannot be processed"
                        + "\n"
                    )
                    exit()
                else:
                    features_dict[composite_key] = element
    all_features = list(set(all_features))
    all_samples = list(set(all_samples))
    return (features_dict, all_features, all_samples)


def table_writer(features_dict, all_features, all_samples, out_name):
    with open(out_name, "w") as output_file:
        header = ["feature_id"] + all_samples
        output_file.write("\t".join(header) + "\n")
        for feature in all_features:
            to_print = []
            to_print.append(feature)
            for sample in all_samples:
                composite_key = (sample, feature)
                if composite_key in features_dict:
                    to_print.append(features_dict[composite_key])
                else:
                    to_print.append("0")
            output_file.write("\t".join(to_print) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="This script integrates multiple count matices into a single output. Each input file is a matrix on tsv format having samples as column names and features as row names in the first column. Provide a list of at least 2 files separated by comma: `matrix_1.tsv matrix_2.tsv matrix_3.tsv`"
    )
    parser.add_argument(
        "--input",
        nargs="+",
        help="List of matrices to be integrated, space separated",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Prefix to be used to name the output file. Default = `integrated_matrix.tsv`",
        required=True,
    )
    args = parser.parse_args()

    if args.output:
        out_name = args.output
    else:
        out_name = "integrated_matrix.tsv"

    ### Calling functions
    features_dict = {}
    all_features = []
    all_samples = []

    for matrix_file in args.input:
        (features_dict, all_features, all_samples) = matrix_parser(
            matrix_file, features_dict, all_features, all_samples
        )

    table_writer(features_dict, all_features, all_samples, out_name)


if __name__ == "__main__":
    main()
