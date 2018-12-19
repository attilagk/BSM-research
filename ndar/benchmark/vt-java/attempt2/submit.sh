#! /usr/bin/env bash

#d=/home/attila/projects/bsm/ndar/benchmark/vt-python/test0
d=`dirname $0`
f1=nichd_btb02.csv
f2=genomics_subject02.csv
f3=genomics_sample03.csv

/usr/lib/jvm/java-8-oracle/bin/java \
    -jar /opt/validationtool/current/validationtool-console-client-1.10.2.jar \
    --validationAPI https://ndar.nih.gov/api/validationtool/v2 \
    --collectionId 2458 \
    --title Benchmark \
    --description "FASTQs and BAMs for Benchmark (CEPH/Utah DNA mixes), Chess lab" \
    --username attilagk \
    --password Chesslab13 \
    --buildPackage \ 
    --submitPackage \
    $f1 $f2 $f3
