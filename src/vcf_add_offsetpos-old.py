#! /usr/bin/env python3

import os
import sys
import vcf
import csv

invcfname = sys.argv[1]
#outvcfname = "offsetpos-" + invcfname
faifname = os.environ["REFSEQ"] + ".fai"

def fai2dict(fai_file):
    """Read a fasta index from path fai_file and create a dictionary with
    contig names as keys and offset as values.
    """
    with open(fai_file, "r") as f:
        csv_reader = csv.reader(f, delimiter = "\t")
        d = {}
        for r in csv_reader:
            d[r[0]] = int(r[2])
        return(d)

def print_records(invcf_file, outvcf_file = invcfname + ".offsetpos", faidict = fai2dict(faifname), nrec = 2, doprint = False):
    """Given an input VCF on path invcf_file and a dictionary faidict
    (produced by the fai2dict function) print the first nrec records' INFO
    field and return the nrec'th record.
    """
    i = 0
    with open(invcf_file, "r") as f, open(outvcf_file, "w") as g:
        invcf = vcf.Reader(f)
        outvcf = vcf.Writer(g, invcf)
        for r in invcf:
            r.INFO["OFFSETPOS"] = r.POS + faidict[r.CHROM]
            outvcf.write_record(r)
            if doprint:
                print(r.INFO)
            if i >= nrec:
                break
            i += 1
        f.close()
        return(r)

if __name__ == "__main__":
    print_records(invcfname, doprint = True).INFO
