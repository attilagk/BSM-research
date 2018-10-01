#! /usr/bin/env python3

import sys
import vcf

invcf = vcf.Reader(open("../results/2018-06-22-variant-meta-caller-test/vmc-prioritize-benchmark/test10/results.vcf", "r"))

#def foo(invcf_file):
#    return(invcf_file * 2)
#    invcf = vcf.Reader(open(invcf_file, "r"))
#    i = 0
#    for record in invcf:
#        i += 1
#        if i >= 10:
#            break
#        print(record.INFO)
#

i = 0
for record in invcf:
    i += 1
    if i >= 10:
        break
    print(record.INFO)


def foo(invcf_file):
    return(invcf_file * 2)
