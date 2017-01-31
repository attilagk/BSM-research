#!/bin/sh

# this script extracts variants that are specific to a given sample

BN=snv-2017-01-24 # basename of the vcf file
#BN=chr1.recode # basename of the vcf file
SAMPLE=Fibro_Common_7 # sample ID

# get positions of all sample-specific variants
vcftools --gzvcf $BN.vcf.gz --keep-filtered PASS --singletons --out $BN

# keep only those that pertain to the given sample
sed "1p; /$SAMPLE/!d" $BN.singletons > $BN.singletons.$SAMPLE

# extract variants given the positions of the sample-specific variants
vcftools --gzvcf $BN.vcf.gz --positions $BN.singletons.$SAMPLE --recode --out $BN.$SAMPLE

# remove header to open file as "tab separated" in Excel
sed '/^##/d' $BN.$SAMPLE.recode.vcf > $BN.$SAMPLE.recode.tsv

# remove header from the original vcf and get the first 1000 variants
gunzip -c $BN.vcf.gz | sed -n '/^#CHROM/,+1000p' > $BN.1-1000.tsv
