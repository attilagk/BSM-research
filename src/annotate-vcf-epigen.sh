#! /usr/bin/env bash

# Annotate VCF with chromatin state from Roadmap Epigenomics project

invcf=$1
outvcf=${2:-"`dirname $invcf`/`basename $invcf .vcf`-chromatin.vcf"}
# Ensure that the following resource path exists
# The name of this BED file ends with "-no-chr-prefix" signifying the fact
# that the chromosome names are without the 'chr' prefix.  It was created from
# the similarly named BED file by a sed substitute operation.
bed=~/projects/bsm/resources/roadmap-epigenomics/chromatin-state/15-state-model/E073_15_coreMarks_stateno-no-chr-prefix.bed
#bed=/big/resources/roadmap-epigenomics/chromatin-state/15-state-model/E073_15_coreMarks_stateno-no-chr-prefix.bed

intersec=`tempfile`
annotations=`tempfile`
headerfile=`tempfile`

# intersect BED and VCF
bedtools intersect -wo -a $bed -b $invcf | cut -f4-6 > $intersec
# create annotations file for bcftools annotate
paste <(cut -f2-3 $intersec) <(cut -f1 $intersec) > $annotations
# bgzipping and indexing
bgzip $annotations
tabix -s1 -b2 -e2 $annotations.gz
# create header
header='##INFO=<ID=ChromatinState,Number=1,Type=Integer,Description="Chromatin state from the Roadmap Epigenomics project">'
echo $header > $headerfile
columns='CHROM,POS,INFO/ChromatinState' 
# annotate VCF
bcftools annotate -a $annotations.gz -c $columns -h $headerfile -o $outvcf $invcf

# clean up
rm $intersec $annotations.gz $annotations.gz.tbi
