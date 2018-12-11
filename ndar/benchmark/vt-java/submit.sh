#! /usr/bin/env bash

d=/home/attila/projects/bsm/ndar/benchmark
f1=$d/nichd_btb02-chess-benchmark-2018-12-10.csv
f2=$d/genomics_subject02-chess-benchmark-2018-12-10.csv
f3=$d/genomics_sample03-chess-benchmark-2018-12-10.csv

file $f1 $f2 $f3; exit

# Oracle 8
# /usr/lib/jvm/java-8-oracle/bin/java \
    # OpenJDK
# java \
/usr/lib/jvm/java-8-oracle/bin/java \
    -jar /opt/validationtool/current/validationtool-console-client-1.10.2.jar \
    --validationAPI https://ndar.nih.gov/api/validationtool/v2 \
    --collectionId 2458 \
    --title chess-benchmark-2018-12-10 \
    --description "DNA Mixes of CEPH/Utah grandparents" \
    --username attilagk \
    --password Chesslab13 \
    --buildPackage \ 
    --submitPackage \
    $f1 $f2 $f3
