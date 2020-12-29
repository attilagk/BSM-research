#! /usr/bin/env bash

#d=/home/attila/projects/bsm/ndar/benchmark/vt-python/test0
d=$1
#d=`realpath $1`
f1=$d/nichd_btb02.csv
f2=$d/genomics_subject02.csv
f3=$d/genomics_sample03.csv

vtcmd \
	$f1 $f2 $f3 \
    -u attilagk \
    -p Chesslab13 \
    -a BSMN-S3 \
    -l $d \
    -t "Chess lab" \
    -d "FASTQs and BAMs of the Chess lab data" -b

exit

    -m $d \
