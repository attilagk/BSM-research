#! /usr/bin/env bash

case $HOSTNAME in
    ada) : ;;
    *) exit ;;
esac

indiv=$1
rownames=$HOME/projects/bsm/tables/summary-rownames.csv

# subject a.k.a. individual
echo $indiv

# BAMs
bamdir=/projects/bsm/alignments/$indiv/
for sample in $(sed -n "/fastq/ { s/.fastq//; p }" $rownames); do
    bam=$bamdir/${indiv}_${sample}.bam
    bai=$bamdir/${indiv}_${sample}.bam.bai
    if test -f $bam && test -f $bai; then
        echo 1
    else
        echo 0
    fi
done
