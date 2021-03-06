#! /usr/bin/env bash

#d=/home/attila/projects/bsm/ndar/benchmark/vt-python/test0
d=`dirname $0`
echo $d; exit
f1=nichd_btb02.csv
f2=genomics_subject02.csv
f3=genomics_sample03.csv

vtcmd \
    --username attilagk \
    --password Chesslab13 \
    --manifestPath $d \
    --collectionID 2458 \
    --title test0 \
    --buildPackage \
    --description "Testing submission" \
    $f1 $f2 $f3
