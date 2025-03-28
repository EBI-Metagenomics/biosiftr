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


def core_parser(core_table):
    all_modules = []
    core_values = {}
    with open(core_table) as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split("\t")
            module_accession = l_line[0]
            completeness = l_line[1]
            pathway_name = l_line[2]
            module_name = module_accession + "|" + pathway_name
            all_modules.append(module_name)
            core_values[module_name] = completeness
    return (all_modules, core_values)


def pan_parser(all_modules, pan_table):
    pan_values = {}
    with open(pan_table) as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split("\t")
            module_accession = l_line[0]
            completeness = l_line[1]
            pathway_name = l_line[2]
            module_name = module_accession + "|" + pathway_name
            all_modules.append(module_name)
            pan_values[module_name] = completeness
    all_modules = list(set(all_modules))
    return (all_modules, pan_values)


def output_writer(all_modules, core_values, pan_values, prefix):
    with open(prefix + "_kegg_comp.tsv", "w") as file_out:
        file_out.write(
            "\t".join(
                [
                    "#module",
                    "pangenome",
                    "core",
                ]
            )
            + "\n"
        )
        for module in all_modules:
            to_print = []
            to_print.append(module)
            if module in pan_values:
                to_print.append(pan_values[module])
            else:
                to_print.append("0")

            if module in core_values:
                to_print.append(core_values[module])
            else:
                to_print.append("0")

            file_out.write("\t".join(to_print) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="This script integrates the output of the kegg completeness tool to build a DB at pangenome level. Provide the prefix to be used in the output file and the core and pangenomic kegg completeness tables"
    )
    parser.add_argument(
        "--core",
        type=str,
        help="Pathways completeness table at core level (*_core.summary.kegg_pathways.tsv)",
        required=True,
    )
    parser.add_argument(
        "--pan",
        type=str,
        help="Pathways completeness table at pangenome level (*_pan.summary.kegg_pathways.tsv)",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Prefix to be used to name the output",
        required=True,
    )
    args = parser.parse_args()

    ### Calling functions
    (all_modules, core_values) = core_parser(args.core)
    (all_modules, pan_values) = pan_parser(all_modules, args.pan)
    output_writer(all_modules, core_values, pan_values, args.output)


if __name__ == "__main__":
    main()
