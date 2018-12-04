#! /usr/bin/env bash

# Given "sequenced-individuals" get sequenced samples and create indiv-sample-fastq-names

# Must be run on Ada.  No command line arguments are required.

case $HOSTNAME in
    ada) : ;;
    *) exit ;;
esac

# store sequenced samples in outfile

outdir=/projects/bsm/attila/results/2018-09-12-sequenced-individuals
outfile=$outdir/sequenced-samples

if test -z $@; then
    set $(sed 's/CMC_\(\S\+\)\t.*$/\1/' $outdir/sequenced-individuals)
fi

for indiv in $@; do
    indivcompact=$(echo $indiv | tr -d '_')
    grep -o "$indiv\s\+\(Control\|SCZ\)" $outdir/sequenced-individuals | tr '\n' '\t'
    find /projects/bsm/reads/ -name "*$indivcompact*.fq.gz" |
        grep -o '\(muscle\|NeuN_pl\|NeuN_mn\)' |
        sort -u | tr '\n' '\t'
    echo
done > $outfile


# for each indiv-sample store fastq file names in indiv-sample-fastq-names

bamdir=/projects/bsm/alignments

cut -f1 $outfile |
    for indiv in `cat`; do
        indivno_=$(echo $indiv | tr -d _ )
        if test -d $bamdir/$indiv; then
            :
        else
            mkdir $bamdir/$indiv
            sed -n "/$indiv/ p" $outfile | cut -f3- |
                for sample in `cat`; do
                    find /projects/bsm/reads/ -name "*${indivno_}_${sample}*.fq.gz" > \
                        $bamdir/$indiv/${indiv}_${sample}-fastq-names
                done
        fi
    done
