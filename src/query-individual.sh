#! /usr/bin/env bash

case $HOSTNAME in
    ada) : ;;
    *) exit ;;
esac

indiv=$1
rownames=$HOME/projects/bsm/tables/summary-rownames.csv

# subject a.k.a. individual
echo $indiv

bamdir=/projects/bsm/alignments/$indiv/
# FASTQs
for sample in $(sed -n "/fastq/ { s/.fastq//; p }" $rownames); do
    fqnames=$bamdir/${indiv}_${sample}-fastq-names
    if test -f $fqnames; then
        echo 1
    else
        echo 0
    fi
done
# BAMs
for sample in $(sed -n "/fastq/ { s/.fastq//; p }" $rownames); do
    bam=$bamdir/${indiv}_${sample}.bam
    bai=$bamdir/${indiv}_${sample}.bam.bai
    if test -f $bam && test -f $bai; then
        echo 1
    else
        echo 0
    fi
done

vcfdir=/projects/bsm/calls/$indiv/snvs
for caller in $(sed -n "/vcf/ { s/.vcf//; p }" $rownames); do
    vcf=$vcfdir/$caller.vcf.gz
    if test -f $vcf; then
        echo 1
    else
        echo 0
    fi
done
