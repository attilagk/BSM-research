#! /usr/bin/env bash

#d=/home/attila/projects/bsm/ndar/benchmark/vt-python/test0
d=/projects/bsm/attila/results/2019-02-19-upload-to-ndar/
f1=$d/CMC_MSSM_106-nichd_btb02_U01MH106891_Chess.csv
f2=$d/CMC_MSSM_106-genomics_subject02_U01MH106891_Chess.csv
f3=$d/CMC_MSSM_106-genomics_sample03_template.csv

vtcmd \
$f1 $f2 $f3 \
-u attilagk \
-p Chesslab13 \
-l '/projects/bsm/reads/', '/projects/bsm/alignments/' \
-a BSMN-S3 \
-t "This is a test" \
-d "Hello, World"
