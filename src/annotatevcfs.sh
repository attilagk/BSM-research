#!/bin/bash

PATH=/home/attila/tools/vcflib/bin:$PATH

SRCDIR="${1:-trial/}"
DESTDIR="${2:-trialres/}"
SCRIPTDIR=$HOME/projects/bsm/src/
REFERENCE="$REFSEQ"

mkdir -p $DESTDIR

for f in $SRCDIR/*.vcf;
do 
    OUT=$DESTDIR/$(basename $f)

    OUT0=$DESTDIR/$(basename $f .vcf).0.vcf
    OUT1=$DESTDIR/$(basename $f .vcf).1.vcf
    OUT2=$DESTDIR/$(basename $f .vcf).2.vcf
    OUT3=$DESTDIR/$(basename $f .vcf).3.vcf

    cp $f $OUT0
    bgzip $OUT0
    tabix -p vcf $OUT0.gz

    echo "#CHROM    POS     REF     ALT     BasesToClosestVariant" > $OUT.vcfdistance.info
    echo "#CHROM    POS     REF     ALT     EntropyLeft_7   EntropyCenter_7 EntropyRight_7" > $OUT.vcfentropy_7.info
    echo "#CHROM    POS     REF     ALT     EntropyLeft_15  EntropyCenter_15        EntropyRight_15" > $OUT.vcfentropy_15.info

    cat $f | vcfdistance | bgzip | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%BasesToClosestVariant\n' >> $OUT.vcfdistance.info
    cat $f | vcfentropy -w 7 -f $REFERENCE | bgzip | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%EntropyLeft\t%EntropyCenter\t%EntropyRight\n' >> $OUT.vcfentropy_7.info
    cat $f | vcfentropy -w 15 -f $REFERENCE | bgzip | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%EntropyLeft\t%EntropyCenter\t%EntropyRight\n' >> $OUT.vcfentropy_15.info

    bgzip $OUT.vcfdistance.info
    tabix -s 1 -b 2 -e 2 -f $OUT.vcfdistance.info.gz

    bgzip $OUT.vcfentropy_7.info
    tabix -s 1 -b 2 -e 2 -f $OUT.vcfentropy_7.info.gz

    bgzip $OUT.vcfentropy_15.info
    tabix -s 1 -b 2 -e 2 -f $OUT.vcfentropy_15.info.gz

    echo "Annotate $f with the number of bases to the closest variant."
    bcftools annotate -a $OUT.vcfdistance.info.gz -c CHROM,POS,REF,ALT,BasesToClosestVariant -h$SCRIPTDIR/vcfdistance.info $OUT0.gz > $OUT1
    bgzip $OUT1
    tabix -p vcf $OUT1.gz

    echo "Annotate $f with the entropy of the reference sequence near the variant in a 7 bp window."
    bcftools annotate -a $OUT.vcfentropy_7.info.gz -c CHROM,POS,REF,ALT,EntropyLeft_7,EntropyCenter_7,EntropyRight_7 -h$SCRIPTDIR/vcfentropy_7.info $OUT1.gz > $OUT2
    bgzip $OUT2
    tabix -p vcf $OUT2.gz

    echo "Annotate $f with the entropy of the reference sequence near the variant in a 15 bp window."
    bcftools annotate -a $OUT.vcfentropy_15.info.gz -c CHROM,POS,REF,ALT,EntropyLeft_15,EntropyCenter_15,EntropyRight_15 -h$SCRIPTDIR/vcfentropy_15.info $OUT2.gz > $OUT3
    bgzip $OUT3
    tabix -p vcf $OUT3.gz

    echo "Left normalize $f"
    bcftools norm -o $OUT -f $REFERENCE -D $OUT3.gz

    rm $OUT0.* $OUT1.* $OUT2.* $OUT3.* $OUT.*

done


