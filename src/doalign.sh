#!/usr/bin/env bash

usage="usage: `basename $0` flowcell lane [library]"
if test $# -lt 2; then
    echo $usage; exit
fi

flowcell=$1
lane=$2
library=${library:-USPD16080279}
indir=~/data/bsm/testdata

#for tissue in muscle NeuN_mn; do
#for tissue in NeuN_pl; do
for tissue in muscle NeuN_mn NeuN_pl; do
    case $tissue in # tissue distinguished by index IDs
        muscle) d=D705;;
        NeuN_pl) d=D706;;
        NeuN_mn) d=D704;;
    esac
    rg="${library}-${d}_${flowcell}_${lane}"
    #rg="ID:${library}-${d}_${flowcell}_${lane} LB:$library SM:MSSM179_$tissue"
    fq1="$indir/MSSM179_${tissue}_${rg}_1*"
    fq2="$indir/MSSM179_${tissue}_${rg}_2*"
    #cmd="align -t6 -rID:$rg -B MSSM179_${tissue}_$rg.bam $fq1 $fq2"
    cmd="align -t6 `parseX10-fnames.1 $fq1` -r PL:Illumina -B MSSM179_${tissue}_$rg.bam $fq1 $fq2"
    $cmd && rm -f $fq1 $fq2
done
