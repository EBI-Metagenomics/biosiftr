#!/usr/bin/env python

import argparse
import pysam
import sys
import re

##### This script process BWA results to compute relative abundance of unique mapped reads
##### Alejandra Escobar, EMBL-EBI
##### Nov 10, 2023


def bam_header( bwa_bam ):
    genomes_len = {}
    bam_file = pysam.AlignmentFile(bwa_bam, "rb")    
    reference_sequences = bam_file.header["SQ"]
    for reference in reference_sequences:
        genome_id = reference["SN"].split('_')[0]
        seq_len = reference["LN"]
        if genome_id in genomes_len:
            added_len = genomes_len[genome_id] + seq_len
            genomes_len[genome_id] = added_len
        else:
            genomes_len[genome_id] = seq_len
    bam_file.close()

    return( genomes_len )


def bam_parser( bam_file ):
    id_thresh = float(90)
    cov_thresh = float(60)
    ani_dicarded = 0
    total_hq = 0
    reads_len_sum = 0
    unique_matches = {}

    with pysam.AlignmentFile(bam_file, 'rb') as input_bam:
        for read in input_bam:
            read_id = str( read.query_name )
            ref_genome = str( read.reference_name ).split('_')[0]
            ani = ( read.query_alignment_length - read.get_tag("NM") ) / float( read.query_alignment_length ) * 100
            cov = read.query_alignment_length / float( read.query_length ) * 100

            # Keeping high-quality mapping reads only
            if all([ ani >= id_thresh, cov >= cov_thresh ]):
                total_hq += 1
                reads_len_sum += read.query_length

                # Unique mapping reads don't have XA tag
                if not 'XA:Z:' in read.tostring():
                    if ref_genome in unique_matches: 
                        unique_matches[ref_genome] += 1
                    else:
                        unique_matches[ref_genome] = 1

            else:
                ani_dicarded += 1

    ave_read_len = float(reads_len_sum)/float(total_hq)

    print('Total number of hq reads: ',total_hq)
    print('Total number of hq unique mapping reads: ',sum(unique_matches.values()))
    print('Total number of low-qual mapping reads: ',ani_dicarded)

    return( unique_matches, ave_read_len )


def FP_control( out_root, genomes_len, unique_matches, ave_read_len ):
    unique_thres01 = []
    cov_threhold = 0.1
    total_unique = 0
    for genome in unique_matches:
        assembly_len = genomes_len[genome]
        mapped_reads = unique_matches[genome]
        genome_coverage = ( float(mapped_reads) * float(ave_read_len) ) / float(assembly_len)
        if genome_coverage >= cov_threhold:
            total_unique = total_unique + mapped_reads
            unique_thres01.append(genome)

    with open(out_root+'.tsv','w') as unique_out:
        unique_out.write('\t'.join([
            'reference_id',
            'reference_len',
            'reads_count',
            'base_cov',
            'rel_abun'
        ])+ '\n' )
        for genome in unique_thres01:
            assembly_len = genomes_len[genome]
            mapped_reads = unique_matches[genome]
            genome_coverage = ( float(mapped_reads) * float(ave_read_len) ) / float(assembly_len)
            relative_abundance = float(mapped_reads) / float(total_unique)
            unique_out.write('\t'.join([
                genome,
                str(assembly_len),
                str(mapped_reads),
                str(genome_coverage),
                str(relative_abundance)
            ])+ '\n' )



def main():
    parser = argparse.ArgumentParser(
        description="This script process the bam file generated using bwa-mem2 with the flag -M. The bam file has to be filtered using -F4 -F256 flags, then sorted and indexed. The output of this script is a table of relative abundance of genomes based on unique mappin greads. Only genomnes having a mean depth coverage >0.1 are considered."
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
        help="To name the output files. Default = u_relab_01",
        required=False,
    )
    args = parser.parse_args()

    if args.prefix:
        out_root = args.prefix
    else:
        out_root = 'u_relab_01'

    ### Calling functions
    ( genomes_len ) = bam_header( args.bwa_bam )

    ( unique_matches, ave_read_len ) = bam_parser( args.bwa_bam ) 

    FP_control( out_root, genomes_len, unique_matches, ave_read_len )


if __name__ == "__main__":
    main()

