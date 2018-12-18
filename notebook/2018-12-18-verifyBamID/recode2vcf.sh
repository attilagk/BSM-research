#! /usr/bin/env bash

dirprefix=/projects/bsm/attila/
wd=$dirprefix/results/2018-12-18-verifyBamID
cd $wd

plink --bfile CMC-preimputed --recode vcf bgz --out CMC-preimputed
