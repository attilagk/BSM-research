#!/usr/bin/env bash

# subsample tumor.bam and run mutect2
usage="usage:$(basename $0) .20 tumor.bam [normal.bam]"

frac=$1
tumorbam=$2
normalbam=${3:-MSSM179_muscle-1MB.bam}
seed=1313
refseq="$HOME/data/GRCh37/karyotypic-order/Homo_sapiens.GRCh37.dna.fa"
# output directory and path to subsampled tumor.bam
outdir="$(basename -s .bam $tumorbam)-$frac"
subsamplebam="$outdir/subsample.bam"
# path to mutect2's output VCF
mutect2dir="$outdir/mutect2"
vcf=$mutect2dir/out.vcf
csv=$mutect2dir/lod.csv

# check syntax
case $frac in
    [10].[0-9][0-9]) : ;;
    *) echo $usage 2>&1; exit 1 ;;
esac

# subsample tumor.bam unless it already exists
if [ ! -e $subsamplebam ]; then
    mkdir $outdir
    case $frac in
        1.00)
            cp $tumorbam $subsamplebam
            samtools index $subsamplebam
            ;;
        *)
            tmpbam=`tempfile -s .bam`
            samtools view -b -s "$seed$frac" $tumorbam > $tmpbam &&
                samtools sort -o $subsamplebam $tmpbam &&
                rm $tmpbam &&
                samtools index $subsamplebam
            ;;
    esac
fi

# run mutect2 on normalbam and subsamplebam unless out.vcf already exists
if [ ! -e $vcf ]; then
    mymutect2 $normalbam $subsamplebam $mutect2dir $refseq
fi

echo "CHROM,POS,NLOD,TLOD" > $csv
sed -n '/^#CHROM/,$ {/^[^#].*/ p}' $vcf |
sed 's/^\([^\t]*\)\t\([^\t]*\)\t.*NLOD=\([.[:digit:]]\+\).*TLOD=\([.[:digit:]]\+\).*$/\1,\2,\3,\4/' >> $csv
