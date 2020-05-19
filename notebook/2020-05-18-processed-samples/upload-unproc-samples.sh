#! /usr/bin/env bash

samplelist=${1:-/projects/bsm/attila/results/2020-05-18-processed-samples/unprocessed_samples}

helper () {
sample=$1
ssample=$(echo $sample | sed 's/_//')
readsdir=/projects/bsm/reads
sampledir=$readsdir/$ssample
echo uploading files from $sampledir
aws s3 cp --recursive $sampledir s3://chesslab-bsmn/Ada/
}

for s in `cat $samplelist`; do
    helper $s
done
