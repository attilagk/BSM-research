#! /usr/bin/env bash

#d=/home/attila/projects/bsm/ndar/benchmark/vt-python/test0
d=`dirname $0`
f1=$d/nichd_btb02.csv
f2=$d/genomics_subject02.csv
f3=$d/genomics_sample03.csv

vtcmd \
	$f1 $f2 $f3 \
    -u attilagk \
    -p Chesslab13 \
    -l /projects/bsm/alignments/ceph-benchmark/ /projects/bsm/reads/2018-01-10-Benchmark-DV-X10/ \
    -a BSMN-S3 \
    -t "U01MH106891, reference_tissue, benchmark/mixin" \
    -d "FASTQs and BAMs for Benchmark (CEPH/Utah DNA mixes), Chess lab" \
    -b
