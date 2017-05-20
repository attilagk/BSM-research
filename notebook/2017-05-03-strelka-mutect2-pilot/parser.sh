#!/usr/bin/env bash

sedcmd='/elapsed/ { s/^\([[:digit:].]\+\)user.*$/\1/; p }' 

echo "seglen,seglen.unit,runtime.strelka,runtime.mutect2"
for l in {1,2,4,8,16}; do
    echo -n $l,${l}MB,
    sed -n "$sedcmd" ${l}MB/muscle-NeuN_pl-strelka.time | tr '\n' ,
    sed -n "$sedcmd" ${l}MB/muscle-NeuN_pl-mutect2.time
done
