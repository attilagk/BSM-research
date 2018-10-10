#! /usr/bin/env python3

import csv
import os
import sys
import vcfpy
import collections

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

def add_offsetpos(invcf_file, faidict = fai2dict(faifname)):
    """Add OFFSETPOS to INFO
    """
    reader = vcfpy.Reader.from_path(invcf_file)
    reader.header.add_info_line(collections.OrderedDict([("ID", "OFFSETPOS"),
        ("Number", 1), ("Type", "Integer"), ("Description", "Position offset by cumulative contig length")]))
    writer = vcfpy.Writer.from_path("/dev/stdout", reader.header)
    for r in reader:
        r.INFO["OFFSETPOS"] = r.POS + faidict[r.CHROM]
        writer.write_record(r)

if __name__ == "__main__":
    add_offsetpos("/dev/stdin")
