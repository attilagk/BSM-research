#! /usr/bin/env python3

import sys
import vcf

invcfname = sys.argv[1]

def print_records(invcf_file, nrec = 2):
    f = open(invcf_file, "r")
    invcf = vcf.Reader(f)
    i = 0
    for r in invcf:
        if i >= nrec:
            break
        print(r.INFO)
        i += 1
    f.close()

if __name__ == "__main__":
    print_records(invcfname)
