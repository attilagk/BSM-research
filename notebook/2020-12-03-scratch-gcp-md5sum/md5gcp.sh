#! /usr/bin/env bash
# get md5 (base64) from Google Cloud Platform store for alignment file

fname=$1
prefix='gs://chesslab-bsmn/alignments/'
uri=$prefix$fname

tmp=`tempfile`

# write gsutil's message into tempfile
# exit with 1 if URI wasn't found
if ! message=$(gsutil hash -m $uri > $tmp); then
	exit 1
fi
# extract md5 from tempfile
cat $tmp | sed -n '/md5/ { s/^.*\s\(\S\+\)\s*$/\1/; p }' && rm $tmp
