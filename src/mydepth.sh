#! /usr/bin/env bash

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
