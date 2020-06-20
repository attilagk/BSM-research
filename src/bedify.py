#! /usr/bin/env python3

# Convert TSV file into BED by setting 1-based start pos (second column) to 0-based.

import sys

def bedify(intsv, outbed):
    with open(intsv) as tsv:
        for line in tsv:
            l = line.split(sep='\t')
            startpos1 = int(l[1])
            startpos0 = startpos1 - 1
            l[1] = startpos0
            print(*l, sep='\t', end='', file=outbed)
    return(None)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('intsv', help='the input TSV file')
    parser.add_argument('-o', '--outbed', help='the output BED file',
            default=sys.stdout)
    parser.add_argument('-l', '--vcflist', help='list of samples and VCF files (/big/results/bsm/calls/filtered-vcfs.tsv)',
            default='/big/results/bsm/calls/filtered-vcfs.tsv')
    args = parser.parse_args()
    bedify(args.intsv, args.outbed)
