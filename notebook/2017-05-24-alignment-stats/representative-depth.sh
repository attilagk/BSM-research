#! /usr/bin/env bash

# Obtain depth at 1 base resolution from a given crhomosome between startpos and endpos

usage="usage: `basename $0` chr [startpos endpos]"
if test $# -eq 0; then
    echo $usage 2>&1; exit 1
fi

chr=${1:-1}
startpos=${2:-50,000,001}
endpos=${3:-52,000,000}

indir=/media/attila/seagate-attila/projects/bsm/MSSM_179
outdir=~/projects/bsm/results/2017-05-24-alignment-stats

for inbam in $indir/*.bam; do
    outbam=$outdir/`basename $inbam`.depth.chr${chr}_$startpos-$endpos
    samtools depth -r $chr:$startpos-$endpos $inbam > $outbam
done
