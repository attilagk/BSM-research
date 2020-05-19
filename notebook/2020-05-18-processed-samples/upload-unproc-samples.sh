#! /usr/bin/env bash
# upload all files to s3://chesslab-bsmn/Ada/ for each sample in sample list

# The assumption is that the files identified by this script are all fq.gz
# files.  This is, however, incorrect because
# /projects/bsmn/reads/MSSM331_NeuN_pl/ contains bam and bam.bai files as
# well.  These were removed by aws s3 rm during the development of
# list-unproc-samples.sh.

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
