#! /usr/bin/env bash

subject=CMC_MSSM_118
d=/projects/bsm/attila/results/2019-02-19-upload-to-ndar/
f1=$d/$subject-nichd_btb02_U01MH106891_Chess.csv
f2=$d/$subject-genomics_subject02_U01MH106891_Chess.csv
f3=$d/$subject-genomics_sample03_U01MH106891_Chess.csv

vtcmd \
$f1 $f2 $f3 \
-u attilagk \
-p Chesslab13 \
-l '/projects/bsm/reads/', '/projects/bsm/alignments/' \
-c 2965 \
-t "This is a test" \
-d "Hello, World"

#-a BSMN-S3 \
