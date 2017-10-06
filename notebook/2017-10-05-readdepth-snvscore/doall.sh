#!/usr/bin/env bash

for f in 1.00 0.50 0.25 0.12 0.06 0.03; do
#for f in 0.06 0.03; do
    ./subsample-bam.sh $f MSSM179_NeuN_pl-1MB.bam MSSM179_muscle-1MB.bam &
done
