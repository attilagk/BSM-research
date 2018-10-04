#! /usr/bin/env python3

import os
import sys
import vcf
import csv

invcfname = sys.argv[1]
if len(sys.argv) <= 2:
    faifname = os.environ["REFSEQ"] + ".fai"
else:
    faifname = sys.argv[2]

def fai2dict(fai_file):
    with open(fai_file, "r") as f:
        csv_reader = csv.reader(f, delimiter = "\t")
        d = {}
        for r in csv_reader:
            d[r[0]] = int(r[2])
        return(d)

def print_records(invcf_file, fai, nrec = 2):
    i = 0
    with open(invcf_file, "r") as f:
        invcf = vcf.Reader(f)
        for r in invcf:
            r.INFO["POS1"] = r.POS + fai[r.CHROM]
            print(r.INFO)
            if i >= nrec:
                break
            i += 1
        f.close()
        return(r)

if __name__ == "__main__":
    print_records(invcfname, fai = fai2dict(faifname))
