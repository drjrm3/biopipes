#!/usr/bin/env python3
""" Utilities for coverage analysis. """

import pysam

#-------------------------------------------------------------------------------
def get_chrom_lens(bamfile, filter_list=[], only_chrom=""):
    """ Get chromosome lengths from bamfile header.
    Args:
        bamfile (pysam.AlignmentFile): Bamfile object from pysam.
        filter_list (list): List of strings to filter out ('alt', 'Un', 'etc').
        only_chrom (str): Only use one chromosome.
    Returns:
        chrom_lens (dict): Dictionary of chromosome names and lengths.
    """

    bam = pysam.AlignmentFile(bamfile, "rb")

    chrom_lens = {}

    # TODO: Better way to go through lengths than converting to str and parsing?
    for line in str(bam.header).splitlines():
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
