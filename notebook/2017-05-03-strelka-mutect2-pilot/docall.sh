#!/usr/bin/env bash

case $1 in
    strelka) caller=strelka;;
    mutect2) caller=mutect2;;
    *) echo "usage: `basename $0` strelka|mutect2 [normal [tumor]]" 2>&1; exit 1;;
esac

normal=${2:-muscle}
tumor=${3:-NeuN_pl}

for seglen in *MB; do
    cd $seglen
    normalbam=MSSM179_$normal-$seglen.bam
    tumorbam=MSSM179_$tumor-$seglen.bam
    >& $normal-$tumor-$caller.time time my$caller $normalbam $tumorbam $normal-$tumor-$caller &
    cd ..
done
