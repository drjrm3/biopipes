#!/usr/bin/env python3
""" Run coverage on bam file.

See https://pysam.readthedocs.io/en/latest/api.html for reference.
"""

import argparse
import pysam
import sys

def get_args():
    """ Get CLI arguments.
    Args: None
    Returns:
        args (argparse.Namespace): Parsed arguments.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-bamfile", type=str, required=True,
                        help="Input bamfile to use.")
    parser.add_argument("-o", "--output-coverage-file", type=str, required=True,
                        help="Output coverage file.")
    parser.add_argument("-w", "--chunksize", type=int, default=100,
                        help="Chunksize to use for parallel computing.")
    parser.add_argument("-p", "--processes", type=int, default=2,
                        help="Number of processes to use.")
    parser.add_argument("-r", "--region", type=str,
                        help="Region to use.")

    args = parser.parse_args()

    return args

def get_chrom_lens(bamfile, filter_list=[]):
    """ Get chromosome lengths from bamfile header.
    Args:
        bamfile (pysam.AlignmentFile): Bamfile object from pysam.
        filter_list (list): List of strings to filter out ('alt', 'Un', 'etc').
    Returns:
        chrom_lens (dict): Dictionary of chromosome names and lengths.
    """

    chrom_lens = {}

    # TODO: Better way to go through lengths than converting to str and parsing?
    for line in str(bamfile.header).splitlines():
        if "@SQ" not in line:
            continue
        words = line.split()
        chrom = words[1].split(':')[1]

        # Filter out strings in filter list from chrom name.
        if any([fstr in chrom for fstr in filter_list]):
            continue
        length = int(words[2].split(':')[1])
        chrom_lens[chrom] = length

    return chrom_lens


def main():
    """ Main entry point. """
    args = get_args()

    bamfile = pysam.AlignmentFile(args.input_bamfile, "rb")
    filter_list=["chrUn", "_alt", "_random"]
    chrom_lens = get_chrom_lens(bamfile, filter_list=filter_list)

    # TODO: Go through 'width' and do 'samtools depth' on each chunksize.
    #       - Is there a pysam module that would work for this?

if __name__ == "__main__":
    main()
