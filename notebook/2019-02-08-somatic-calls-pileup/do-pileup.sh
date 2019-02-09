#! /usr/bin/env bash

# To be run on Ada

incsv=/projects/bsm/attila/results/2018-04-08-call-set-concordance/commonsample-isec.csv
wd=/projects/bsm/attila/results/2019-02-08-somatic-calls-pileup
out=$wd/results.pileup
poslist=$wd/position-list
bamdir=/projects/bsm/alignments/MSSM_106
bams=$bamdir/MSSM_106_*.bam

echo "bam	tissue	chromosome	position	reference	depth	bases base_qual" > $out

for bam in $bams; do
    tissue=$(echo $bam | sed 's/^.*MSSM_106_\(NeuN_mn\|NeuN_pl\|muscle\)\.bam/\1/')
    for region in X:140336646-140336646 X:1355567-1355567; do
        pileup=$(samtools mpileup -f $REFSEQ -r $region $bam)
        echo -e "$bam\t$tissue\t$pileup" >> $out
    done
done
