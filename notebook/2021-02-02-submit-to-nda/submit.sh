#! /usr/bin/env bash
# make asynchronous submissions

cd $HOME/projects/bsm/results/2021-02-02-submit-to-nda
validate="vtcmd -b -t title -d description -s3 nda-bsmn-scratch -u $NDA_USER -p $NDA_PASSWORD -c 2965 -w"
manifests="nichd_btb02-2021-02-02.csv genomics_subject02-2021-02-02.csv"
dftypes="flagstat.txt ploidy_2.vcf.gz ploidy_2.vcf.gz.tbi ploidy_12.vcf.gz ploidy_12.vcf.gz.tbi ploidy_50.vcf.gz ploidy_50.vcf.gz.tbi unmapped.bam"

for prefix in 3338602 3340241; do
	gsublist="genomics_sample03-${prefix}*.csv"
	for gsub in $gsublist; do
		echo '===='
		echo $gsub
		$validate -pre $prefix $gsub $manifests 1> ${prefix}-${gsub}.log 2>&1 &
		echo '----'
	done
done
