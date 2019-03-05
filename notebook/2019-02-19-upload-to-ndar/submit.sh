#! /usr/bin/env bash

subject=$1
Dx=$2
dobuild=${3:-nobuild}

if test $dobuild == '-b'; then
    buildstr=$dobuild
elif test $dobuild == nobuild; then
    buildstr=''
else
    echo wrong argument; exit
fi

grant=U01MH106891
collection=2965
title="$Dx subject $subject"
description="Unmapped and mapped reads from bulk sequencing for $Dx subject $subject; NIH $grant; Chess Lab, Mount Sinai, New York"

d=/projects/bsm/attila/results/2019-02-19-upload-to-ndar/
f1=$d/$subject-nichd_btb02_${grant}_Chess.csv
f2=$d/$subject-genomics_subject02_${grant}_Chess.csv
f3=$d/$subject-genomics_sample03_${grant}_Chess.csv

vtcmd \
$f1 $f2 $f3 \
-u andrewchess \
-p Bern1e2017 \
-l '/projects/bsm/' \
-c $collection \
-t "$title" \
-d "$description" $buildstr
