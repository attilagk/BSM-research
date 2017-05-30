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
    paste $tmp1 $tmp2 > $indir/callset-sizes.tsv && rm $tmp1 $tmp2
}

maindir=$HOME/projects/bsm/results/2017-05-29-vcf-comparisons
for t in snvs indels; do
    dosummary $maindir/2_cmp-reftissues/$t/
    for ref in NeuN_mn muscle; do
        dosummary $maindir/1_isec-callers/$ref-NeuN_pl/$t/
    done
done
