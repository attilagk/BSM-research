#! /usr/bin/env bash

# create sample_list.txt of unprocessed samples for bsmn-pipeline

resdir=/projects/bsm/attila/results/2020-05-18-processed-samples
slist=${1:-$resdir/unprocessed_samples}
sample_list=$resdir/sample_list.txt # output is the input of bsmn-pipeline


helper () {
sample=$1
ssample=$(echo $sample | sed 's/_//')
readsdir=/projects/bsm/reads
sampledir=$readsdir/$ssample
for fq in $sampledir/*.fq.gz; do
    bn=$(basename $fq)
    echo -e "$sample\t$bn\ts3://chesslab-bsmn/Ada/$bn"
done
}

echo -e "#sample_id\tfile_name\tlocation" > $sample_list
for s in `cat $slist`; do
    helper $s
done >> $sample_list
