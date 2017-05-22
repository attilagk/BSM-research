#!/usr/bin/env bash

sedcmd='/elapsed/ { s/^\([[:digit:].]\+\)user.*$/\1/; p }' 

echo "seglen,seglen.unit,caller,normal,tumor,tissue.pair,real.time"
#for l in 1 2 4; do
for l in {1,2,4,8,16,32}; do
    for caller in strelka mutect2; do
        for f in ${l}MB/*-$caller.time; do
            echo -n $(($l * 1000000)),${l}MB,$caller,
            sed 's/^\([^-]\+\)-\([^-]\+\)-.*$/\1,\2,\1-\2/' <<< `basename $f` | tr '\n' ,
            sed -n "$sedcmd" $f
        done
    done
done
