#!/usr/bin/env python

import argparse
import pysam
import sys
import re

##### This script process BWA results to compute per genome coverage
##### Alejandra Escobar, EMBL-EBI
##### v1.0 Nov 10, 2023
##### v1.1 Mar 13, 2024
#####   Function to calculate relative abundance added


def bam_header(bwa_bam):
    genomes_len = {}
    bam_file = pysam.AlignmentFile(bwa_bam, "rb")
    reference_sequences = bam_file.header["SQ"]
    for reference in reference_sequences:
        genome_id = reference["SN"].split("_")[0]
        seq_len = reference["LN"]
        if genome_id in genomes_len:
            added_len = genomes_len[genome_id] + seq_len
            genomes_len[genome_id] = added_len
        else:
            genomes_len[genome_id] = seq_len
    bam_file.close()

    return genomes_len


def bam_parser(bam_file):
    id_thresh = float(90)
    cov_thresh = float(60)
    ani_dicarded = 0
    total_hq = 0
    reads_len_sum = 0
    multimapping_counter = 0
    all_matches = {}
    unique_matches = {}
    best_matches = {}
    multi_matches = {}

    with pysam.AlignmentFile(bam_file, "rb") as input_bam:
        for read in input_bam:
            read_id = str(read.query_name)
            ref_genome = str(read.reference_name).split("_")[0]
            ani = (
                (read.query_alignment_length - read.get_tag("NM"))
                / float(read.query_alignment_length)
                * 100
            )
            cov = read.query_alignment_length / float(read.query_length) * 100

            # Keeping high-quality mapping reads only
            if all([ani >= id_thresh, cov >= cov_thresh]):
                total_hq += 1
                reads_len_sum += read.query_length

                if ref_genome in all_matches:
                    all_matches[ref_genome] += 1
                else:
                    all_matches[ref_genome] = 1

                # Catching multimapping reads
                if "XA:Z:" in read.tostring():
                    multimapping_counter += 1
                    if ref_genome in best_matches:
                        best_matches[ref_genome] += 1
                    else:
                        best_matches[ref_genome] = 1

                    if ref_genome in multi_matches:
                        multi_matches[ref_genome] += 1
                    else:
                        multi_matches[ref_genome] = 1

                    extra_matches = read.get_tag("XA")
                    for alignment in extra_matches.split(";"):
                        if len(alignment) > 0:
                            target, start, cigar, mapq = alignment.split(",")
                            matches = sum(
                                int(length)
                                for length, op in re.findall(r"(\d+)([M=X])", cigar)
                            )
                            xa_cov = float(matches) / float(read.query_length) * 100
                            xa_genome = target.split("_")[0]
                            if xa_cov > 60:
                                if xa_genome in multi_matches:
                                    multi_matches[xa_genome] += 1
                                else:
                                    multi_matches[xa_genome] = 1

                                # Including multi-mapping reads on all matches counts
                                if ref_genome in all_matches:
                                    all_matches[ref_genome] += 1
                                else:
                                    all_matches[ref_genome] = 1

                else:
                    # Unique mapping reads don't have XA tag
                    if ref_genome in unique_matches:
                        unique_matches[ref_genome] += 1
                    else:
                        unique_matches[ref_genome] = 1

                    # Including unique matches on best matches counter
                    if ref_genome in best_matches:
                        best_matches[ref_genome] += 1
                    else:
                        best_matches[ref_genome] = 1

            else:
                ani_dicarded += 1

    ave_read_len = float(reads_len_sum) / float(total_hq)
    print("Total number of primary alignments discarded by ani: ", ani_dicarded)

    return (
        all_matches,
        unique_matches,
        best_matches,
        multi_matches,
        ave_read_len,
        multimapping_counter,
    )


def coverage_calc(
    out_root,
    genomes_len,
    all_matches,
    unique_matches,
    best_matches,
    multi_matches,
    ave_read_len,
    multimapping_counter,
):
    # Computing coverage per genome for unique alignments
    with open(out_root + "_unique.tsv", "w") as unique_out:
        unique_out.write(
            "\t".join(
                ["reference_id", "reference_len", "reads_count", "base_cov", "rel_abun"]
            )
            + "\n"
        )
        total_unique_reads = total_reads_counter(unique_matches)
        for genome in unique_matches:
            assembly_len = genomes_len[genome]
            mapped_reads = unique_matches[genome]
            genome_coverage = (float(mapped_reads) * float(ave_read_len)) / float(
                assembly_len
            )
            relative_abundance = float(mapped_reads) / float(total_unique_reads)
            unique_out.write(
                "\t".join(
                    [
                        genome,
                        str(assembly_len),
                        str(mapped_reads),
                        str(genome_coverage),
                        str(relative_abundance),
                    ]
                )
                + "\n"
            )

    # Computing coverage per genome for unique plus best alignments
    with open(out_root + "_best.tsv", "w") as best_out:
        best_out.write(
            "\t".join(
                ["reference_id", "reference_len", "reads_count", "base_cov", "rel_abun"]
            )
            + "\n"
        )
        total_best_reads = total_reads_counter(best_matches)
        for genome in best_matches:
            assembly_len = genomes_len[genome]
            mapped_reads = best_matches[genome]
            genome_coverage = (float(mapped_reads) * float(ave_read_len)) / float(
                assembly_len
            )
            relative_abundance = float(mapped_reads) / float(total_best_reads)
            best_out.write(
                "\t".join(
                    [
                        genome,
                        str(assembly_len),
                        str(mapped_reads),
                        str(genome_coverage),
                        str(relative_abundance),
                    ]
                )
                + "\n"
            )

    # Computing coverage per genome for all alignments
    with open(out_root + "_all.tsv", "w") as all_out:
        all_out.write(
            "\t".join(
                ["reference_id", "reference_len", "reads_count", "base_cov", "rel_abun"]
            )
            + "\n"
        )
        total_all_reads = total_reads_counter(all_matches)
        for genome in all_matches:
            assembly_len = genomes_len[genome]
            mapped_reads = all_matches[genome]
            genome_coverage = (float(mapped_reads) * float(ave_read_len)) / float(
                assembly_len
            )
            relative_abundance = float(mapped_reads) / float(total_all_reads)
            all_out.write(
                "\t".join(
                    [
                        genome,
                        str(assembly_len),
                        str(mapped_reads),
                        str(genome_coverage),
                        str(relative_abundance),
                    ]
                )
                + "\n"
            )

    # Reporting multimapping stats
    with open(out_root + "_multimap_report.tsv", "w") as report_out:
        targets_num = len(list(multi_matches.keys()))
        report_out.write(
            "Total number of targets with multimapping reads: "
            + str(targets_num)
            + "\n"
        )
        report_out.write(
            "Total number of multimapping reads: " + str(multimapping_counter) + "\n"
        )

        sorted_dict = dict(
            sorted(multi_matches.items(), key=lambda item: item[1], reverse=True)
        )
        for genome in sorted_dict:
            report_out.write(genome + "\t" + str(sorted_dict[genome]) + "\n")


def total_reads_counter(matches_dir):
    total_mapped_reads = 0
    for genome in matches_dir:
        mapped_reads = matches_dir[genome]
        total_mapped_reads = total_mapped_reads + mapped_reads
    return total_mapped_reads


def main():
    parser = argparse.ArgumentParser(
        description="This script process the bam file generated using bwa-mem2 with the flag -M. The bam file has to be filtered using -F4 -F256 flags, then sorted and indexed. The outputs of this script are three tables with the base (mean depth) coverage per genome using all mapped reads, unique mapping reads, unique + best mapped reads. Stats on the multimaping targets are generated as well."
    )
    parser.add_argument(
        "--bwa_bam",
        type=str,
        help="bwa-mem2 output in bam format filtered, sorted, and indexed",
        required=True,
    )
    parser.add_argument(
        "--prefix",
        type=str,
        help="To name the output files. Default = coverage",
        required=False,
    )
    args = parser.parse_args()

    if args.prefix:
        out_root = args.prefix
    else:
        out_root = "coverage"

    ### Calling functions
    (genomes_len) = bam_header(args.bwa_bam)

    (
        all_matches,
        unique_matches,
        best_matches,
        multi_matches,
        ave_read_len,
        multimapping_counter,
    ) = bam_parser(args.bwa_bam)

    coverage_calc(
        out_root,
        genomes_len,
        all_matches,
        unique_matches,
        best_matches,
        multi_matches,
        ave_read_len,
        multimapping_counter,
    )


if __name__ == "__main__":
    main()
