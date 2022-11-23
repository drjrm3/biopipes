#!/usr/bin/env python3

import argparse

def get_args():
    """ Get CLI args. """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-chrlens", type=str, required=True)
    parser.add_argument("-o", "--output-bed", type=str, required=True)
    parser.add_argument("-b", "--bin", type=int, default=100)

    args = parser.parse_args()

    return args

def main():
    """ Main routine. """
    args = get_args()

    with open(args.output_bed, "w") as fout:
        with open(args.input_chrlens, "r") as finp:
            line = finp.readline().strip()
            while line:
                words = line.split()
                chrom = words[0]
                length = int(words[1])
                base = 0
                while base < length + args.bin:
                    start = base
                    end = base + args.bin - 1
                    print("%s\t%i\t%i"%(chrom, start, end), file=fout)
                    base += args.bin
                line = finp.readline().strip()


if __name__ == "__main__":
    main()

