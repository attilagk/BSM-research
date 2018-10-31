#!/usr/bin/env bash

# Annotate using vcflib a single VCF.
# Input VCF from STDIN, output VCF to STDOUT in uncompressed (default) or
# compressed output type.

usage="`realpath $0` < invcf > outvcf"

PATH=/home/attila/tools/vcflib/bin:$PATH
outputType=v

while getopts "O:" opt; do
    case $opt in
        O) outputType=$OPTARG;;
    esac
done
shift $(($OPTIND - 1))

SCRIPTDIR=$HOME/projects/bsm/src/
REFERENCE="$REFSEQ"


invcf=`tempfile`
bcftools view -O v -o $invcf
outvcf=`tempfile`
OUT0=`tempfile`
OUT1=`tempfile`
OUT2=`tempfile`
OUT3=`tempfile`

bcftools view -O v -o $OUT0 $invcf
bgzip $OUT0
tabix -p vcf $OUT0.gz

echo "#CHROM    POS     REF     ALT     BasesToClosestVariant" > $outvcf.vcfdistance.info
echo "#CHROM    POS     REF     ALT     EntropyLeft_7   EntropyCenter_7 EntropyRight_7" > $outvcf.vcfentropy_7.info
echo "#CHROM    POS     REF     ALT     EntropyLeft_15  EntropyCenter_15        EntropyRight_15" > $outvcf.vcfentropy_15.info

cat $invcf | vcfdistance | bgzip | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%BasesToClosestVariant\n' >> $outvcf.vcfdistance.info
cat $invcf | vcfentropy -w 7 -f $REFERENCE | bgzip | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%EntropyLeft\t%EntropyCenter\t%EntropyRight\n' >> $outvcf.vcfentropy_7.info
cat $invcf | vcfentropy -w 15 -f $REFERENCE | bgzip | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%EntropyLeft\t%EntropyCenter\t%EntropyRight\n' >> $outvcf.vcfentropy_15.info

bgzip $outvcf.vcfdistance.info
tabix -s 1 -b 2 -e 2 -f $outvcf.vcfdistance.info.gz

bgzip $outvcf.vcfentropy_7.info
tabix -s 1 -b 2 -e 2 -f $outvcf.vcfentropy_7.info.gz

bgzip $outvcf.vcfentropy_15.info
tabix -s 1 -b 2 -e 2 -f $outvcf.vcfentropy_15.info.gz

1>&2 echo "Annotate $invcf with the number of bases to the closest variant."
bcftools annotate -a $outvcf.vcfdistance.info.gz -c CHROM,POS,REF,ALT,BasesToClosestVariant -h$SCRIPTDIR/vcfdistance.info $OUT0.gz > $OUT1
bgzip $OUT1
tabix -p vcf $OUT1.gz

1>&2 echo "Annotate $invcf with the entropy of the reference sequence near the variant in a 7 bp window."
bcftools annotate -a $outvcf.vcfentropy_7.info.gz -c CHROM,POS,REF,ALT,EntropyLeft_7,EntropyCenter_7,EntropyRight_7 -h$SCRIPTDIR/vcfentropy_7.info $OUT1.gz > $OUT2
bgzip $OUT2
tabix -p vcf $OUT2.gz

1>&2 echo "Annotate $invcf with the entropy of the reference sequence near the variant in a 15 bp window."
bcftools annotate -a $outvcf.vcfentropy_15.info.gz -c CHROM,POS,REF,ALT,EntropyLeft_15,EntropyCenter_15,EntropyRight_15 -h$SCRIPTDIR/vcfentropy_15.info $OUT2.gz > $OUT3
bgzip $OUT3
tabix -p vcf $OUT3.gz

1>&2 echo "Left normalize $invcf"
bcftools norm -O $outputType -f $REFERENCE -D $OUT3.gz

#rm $OUT0.* $OUT1.* $OUT2.* $OUT3.* $outvcf.* $invcf
