#! /usr/bin/env bash

vcf=$1
outvcf=${2:-"`dirname $vcf`/`basename $vcf .vcf`-epigen.vcf"}
bed=/big/resources/roadmap-epigenomics/chromatin-state/15-state-model/E073_15_coreMarks_stateno-no-chr-prefix.bed

intersec=`tempfile`
annotations=`tempfile`
headerfile=`tempfile`

# intersect BED and VCF
bedtools intersect -wo -a $bed -b $vcf | cut -f4-6 > $intersec
# create annotations file for bcftools annotate
paste <(cut -f2-3 $intersec) <(cut -f1 $intersec) > $annotations
# bgzipping and indexing
bgzip $annotations
tabix -s1 -b2 -e2 $annotations.gz
# create header
header='##INFO=<ID=EPIGEN_STATE,Number=1,Type=Integer,Description="Epigenetic state from the Roadmap Epigenomics project">'
echo $header > $headerfile
columns='CHROM,POS,INFO/EPIGEN_STATE' 
# annotate VCF
bcftools annotate -a $annotations.gz -c $columns -h $headerfile -o $outvcf $vcf

# clean up
rm $intersec $annotations.gz $annotations.gz.tbi
