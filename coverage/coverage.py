#!/usr/bin/env python3
""" Run coverage on bam file.

See https://pysam.readthedocs.io/en/latest/api.html for reference.
"""

import argparse
import pysam
import multiprocessing
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
    parser.add_argument("-c", "--chrom", type=str, required=True,
                        help="Chromosome to calculate on.")
    parser.add_argument("-p", "--processes", type=int, default=2,
                        help="Number of processes to use.")

    args = parser.parse_args()

    return args

def get_chrom_lens(bamfile, filter_list=[], only_chrom=""):
    """ Get chromosome lengths from bamfile header.
    Args:
        bamfile (pysam.AlignmentFile): Bamfile object from pysam.
        filter_list (list): List of strings to filter out ('alt', 'Un', 'etc').
        only_chrom (str): Only use one chromosome.
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
        if only_chrom and only_chrom != chrom:
            continue

        # Filter out strings in filter list from chrom name.
        if any([fstr in chrom for fstr in filter_list]):
            continue
        length = int(words[2].split(':')[1])
        chrom_lens[chrom] = length

    return chrom_lens

def samtools_depth(args):
    """ Call samtools depth {args}.
    Args:
        args (str): Arguments to pass to samtools depth.
    Returns:
        returncode (int): Returncode of call.
        outlines (str): Stdout.
        errlines (str): Stderr.
    """

    return 0, "OUT: " + args, "ERR: " + args

def gen_coverage(bamfile, outfile, chrom_lens, processes):
    """ Generate coverage file *using samtools*. Note - move to using pysam
        if it helps speed or is easier.
    Args:
        bamfile (str): Input bamfile to use.
        outfile (str): Output coverage file to write to.
        processes (int): Number of processes to use.
        chrom_lens (dict): Dictionary of {chrom: length}.
    Returns: None
    """

    # For now set all regions to just chromosomes.
    regions = chrom_lens.keys()

    args = []
    for region in ["chrM", "chr9"]:
        args.append([f"-r {region} {bamfile}"])

    with multiprocessing.Pool(processes=processes) as pool:
        output = pool.starmap(samtools_depth, args)
        print(output)









def main():
    """ Main entry point. """
    args = get_args()

    print(pysam.__file__)

    bamfile = pysam.AlignmentFile(args.input_bamfile, "rb")
    filter_list=["chrUn", "_alt", "_random"]
    # TODO: Pass args.chrom to get_chrom_lens and only generate {chrom: length}.
    chrom_lens = get_chrom_lens(bamfile, filter_list=filter_list)

    gen_coverage(args.input_bamfile, args.output_coverage_file, chrom_lens,
                 args.processes)

    # TODO: Go through 'width' and do 'samtools depth' on each chunksize.
    #       - Is there a pysam module that would work for this?

if __name__ == "__main__":
    main()
