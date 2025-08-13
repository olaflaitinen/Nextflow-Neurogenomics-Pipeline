#!/usr/bin/env python

import os
import sys
import argparse

def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample,fastq_1,fastq_2,status
    """
    sample_mapping_dict = {}
    with open(file_in, "r") as fin:
        ## Check header
        MIN_COLS = 4
        HEADER = ["sample", "fastq_1", "fastq_2", "status"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[:len(HEADER)] != HEADER:
            print(f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}")
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns
            if len(lspl) < len(HEADER):
                print_error(f"Invalid number of columns (minimum = {len(HEADER)})!", "Line", line)

            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(f"Invalid number of populated columns (minimum = {MIN_COLS})!", "Line", line)

            ## Check sample name entries
            sample, fastq_1, fastq_2, status = lspl[:MIN_COLS]
            if sample.find(" ") != -1:
                print(f"WARNING: Sample name contains spaces. Splicing to '{sample.split()[0]}'.")
                sample = sample.split()[0]
            if not sample:
                print_error("Sample name is missing!", "Line", line)

            ## Auto-detect paired-end/single-end
            if sample and fastq_1 and fastq_2:  # Paired-end
                sample_info = [sample, fastq_1, fastq_2, status]
            else:
                print_error("Invalid FastQ setup. Requires paired-end reads.", "Line", line)

            ## Create sample mapping dictionary
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                print_error("Samplesheet contains duplicate sample name!", "Line", line)

    ## Write validated samplesheet to file
    with open(file_out, "w") as fout:
        fout.write(",".join(["sample", "fastq_1", "fastq_2", "status"]) + "\n")
        for sample in sorted(sample_mapping_dict.keys()):
            for val in sample_mapping_dict[sample]:
                fout.write(",".join(val) + "\n")

def print_error(error, context, line=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if line:
        error_str += f"\n       {context.strip()}: '{line.strip()}'"
    sys.exit(error_str)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Sanity check and format a samplesheet for the pipeline.""")
    parser.add_argument("FILE_IN", help="Input samplesheet.")
    parser.add-argument("FILE_OUT", help="Output file.")
    args = parser.parse_args()
    check_samplesheet(args.FILE_IN, args.FILE_OUT)
