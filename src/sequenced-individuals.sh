#! /usr/bin/env bash

# Lists CMC individuals for which there are FASTQ files

indivdir=/projects/bsm/data/dnalib
indivs=$indivdir/cmc-individuals
readsdir=/projects/bsm/reads
for indiv in `cat $indivs`; do
    indiv="$(echo $indiv | sed 's/CMC_//' | tr -d '_')"
    if test $(find $readsdir -name "${indiv}*.fq*" | wc -w) -ne 0; then
        echo $indiv
    fi
done
