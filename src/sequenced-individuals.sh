#! /usr/bin/env bash

# Lists CMC individuals for which there are FASTQ files

indivdir=/projects/bsm/data/dnalib
indivs=$indivdir/cmc-individuals
chesscsv=$indivdir/BSM_Project_Chess.csv
sortedunique=`mktemp`
cut -f1,12 -d, $chesscsv | sort -u > $sortedunique
readsdir=/projects/bsm/reads
for indiv in `cat $indivs`; do
    indivshort="$(echo $indiv | sed 's/CMC_//' | tr -d '_')"
    if test $(find $readsdir -name "${indivshort}*.fq*" | wc -w) -ne 0; then
        #echo $indiv
        grep $indiv $sortedunique | tr ',' '\t'
    fi
done
rm $sortedunique
