#!/usr/bin/env python3
""" Run coverage on bam file.

See https://pysam.readthedocs.io/en/latest/api.html for reference.
"""

import argparse
import multiprocessing
import shlex
import subprocess

from coverage_utils import get_chrom_lens

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
    parser.add_argument("-c", "--chrom", type=str,
                        help="Chromosome to calculate on.")
    parser.add_argument("-p", "--processes", type=int, default=2,
                        help="Number of processes to use.")

    args = parser.parse_args()

    return args


def samtools_depth(args):
    """ Call samtools depth {args}.
    Args:
        args (str): Arguments to pass to samtools depth.
    Returns:
        returncode (int): Returncode of call.
        outlines (str): Stdout.
        errlines (str): Stderr.
    """

    command = shlex.split(f"samtools depth {args}")

    proc = subprocess.run(command, shell=False, check=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

    stdout_str = ""
    if proc.stdout:
        stdout_str = proc.stdout.decode("utf-8")
    return_code = proc.returncode

    return return_code, stdout_str

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

    # Build samtools depth commands list.
    args = []
    for region in regions:
        args.append([f"-@ 2 -r {region} {bamfile}"])

    # Run them all.
    with multiprocessing.Pool(processes=processes) as pool:
        procs = pool.starmap(samtools_depth, args)
        with open(outfile, "w") as fout:
            for proc in procs:
                print(proc[1], file=fout)


def main():
    """ Main entry point. """
    args = get_args()

    #bamfile = pysam.AlignmentFile(args.input_bamfile, "rb")
    filter_list=["chrUn", "_alt", "_random"]

    # TODO: Pass args.chrom to get_chrom_lens and only generate {chrom: length}.
    chrom_lens = get_chrom_lens(args.input_bamfile, filter_list=filter_list)

    gen_coverage(args.input_bamfile, args.output_coverage_file, chrom_lens,
                 args.processes)


if __name__ == "__main__":
    main()
