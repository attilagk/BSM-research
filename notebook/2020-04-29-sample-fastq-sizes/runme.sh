#! /usr/bin/env bash

rundir=$(dirname $(realpath $0))
$rundir/get-fastq-nblocks
find /projects/bsm/alignments/{MSSM,PITT}* -name '*bam' | xargs $rundir/get_nreads.py
