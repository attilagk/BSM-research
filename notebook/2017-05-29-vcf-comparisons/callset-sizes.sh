#! /usr/bin/env bash

dosummary () {
    indir=$1
    tmp1=`mktemp`
    tmp2=`mktemp`
    for v in 000{0..2}.vcf; do
        grep -v '^##' $indir/$v | wc -l
        readme=$indir/README.txt
    done > $tmp1
    sed -e \
    "1,/^Using/ d;
    s|$indir||g;
    s/for records.*\(private to.*$\|shared by.*$\)/\1/;
    /0003\.vcf/ d" $readme > $tmp2
    paste $tmp1 $tmp2 > $indir/summary.tsv && rm $tmp1 $tmp2
}

maindir=$HOME/projects/bsm/results/2017-05-29-vcf-comparisons
dosummary $maindir/cmp-reftissues/indels/
dosummary $maindir/caller-isec/muscle-NeuN_pl/indels/
dosummary $maindir/cmp-reftissues/snps/
dosummary $maindir/caller-isec/muscle-NeuN_pl/snps/
