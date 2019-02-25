#! /usr/bin/env bash

# Summarize progress in workflow on CMC individuals

usage="$0"

# only run on Ada
case $HOSTNAME in
    ada) : ;;
    *) exit ;;
esac

# assignments
tablesdir=$HOME/projects/bsm/tables
rownames=$tablesdir/summary-rownames.csv
resultdir=/projects/bsm/attila/results/2018-09-12-sequenced-individuals
summarydir=$resultdir/summary-cmc
summaryfile=$resultdir/summary-cmc.tsv
if ! test -d $summarydir; then
    mkdir $summarydir
fi
seqind=/projects/bsm/attila/results/2018-09-12-sequenced-individuals/sequenced-individuals

# get info for a single indiv
do1indiv () {
    # subject a.k.a. individual
    indiv=$1
    bamdir=/projects/bsm/alignments/$indiv/
    # subject
    echo $indiv
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
    # VCFs
    for spair in $samplepairs; do
        vcfdir=/projects/bsm/calls/$indiv/$spair/snvs
        for caller in $(sed -n "/^mu.*vcf/ { s/^\(mu\|ne\)_\(.\+\).vcf/\2/; p }" $rownames); do
            vcf=$vcfdir/$caller.vcf.gz
            if test -f $vcf; then
                echo 1
            else
                echo 0
            fi
        done
    done
}

# do all individuals
cp $rownames $summaryfile
indivs="$(cut -f1 $seqind | sed 's/CMC_//')"
samplepairs="NeuN_pl-muscle NeuN_pl-NeuN_mn"
for ind in $indivs; do
    do1indiv $ind > $summarydir/$ind
done
cd $summarydir
paste $rownames $indivs > $summaryfile
cd .. && rm -r $summarydir
