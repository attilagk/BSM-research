#! /usr/bin/env bash

# Get depth subsampled at every period'th (default: every 1000th) nucleotide.
usage="`basename $0` [-p 1000] [-t targetdir ] aln1.bam aln2.bam ..."

# Depth is written into targetdir with filenames such as aln1.bam.depth.1000

period=1000 # default
while getopts ":p:t:" opt; do
    case $opt in
        p) period=$OPTARG;;
        t) targetdir=$OPTARG;;
    esac
done
shift $(( $OPTIND - 1))

for bampath; do
    dest=`dirname $bampath`
    bam=`basename $bampath`
    tdir=${targetdir:-$dest}
    samtools depth -a $bampath | sed -n "1~$period p" > "$tdir/$bam.depth.$period"
done
