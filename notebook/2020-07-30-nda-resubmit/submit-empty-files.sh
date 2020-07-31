#! /usr/bin/env bash

csv=`realpath $1`
cd `dirname $csv`

sed -n '2,$ {s/,.*$//; p}' $csv | for fpath in `cat`; do
	dirn=$(echo $fpath | sed 's/\/[^/]\+$//')
	if ! test -e $fpath; then
		# take care of directory
		if ! test -z $dirn; then
			# make directory unless it exists
			if ! test -d $dirn; then mkdir -p $dirn; fi
		fi
		# create empty file
		touch $fpath;
	fi
done

vtcmd nichd_btb02.csv genomics_subject02.csv genomics_sample03.csv \
    -u attilagk \
    -p Chesslab13 \
    -a BSMN-S3 \
    -t "Chess lab data" \
    -d "FASTQs and BAMs for the BSMN project, Chess lab, Mount Sinai" \
    -b
