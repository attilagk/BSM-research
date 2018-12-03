#! /usr/bin/env bash

case $HOSTNAME in
    ada) : ;;
    *) exit ;;
esac

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
