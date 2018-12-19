#! /usr/bin/env bash

#d=/home/attila/projects/bsm/ndar/benchmark/vt-python/test0
d=`dirname $0`
f1=nichd_btb02.csv
f2=genomics_subject02.csv
f3=genomics_sample03.csv

vtcmd \
	$f1 $f2 $f3 \
    -u attilagk \
    -p Chesslab13 \
    -l /projects/bsm /projects/bsm/reads/2018-01-10-Benchmark-DV-X10 \
    -c 2458 \
    -t Benchmark \
    -d "FASTQs and BAMs for Benchmark (CEPH/Utah DNA mixes), Chess lab" \
    -b
