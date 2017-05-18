#!/usr/bin/env bash

getopts "d:" opt
destdir=$OPTARG
shift $(($OPTIND - 1))
starttime=`date -Iseconds`

fix1bam () {
    inbam=$1
    path_bn="$destdir/`basename $inbam .bam`"
    rgfix0="${path_bn}-rgfix0.bam"
    rgfix1="${path_bn}-rgfix1.bam"
    #rgfix0=`mktemp`
    #rgfix1=`mktemp`
    outbam="${path_bn}.bam"
    # fix read groups with samtools
    samtools addreplacerg `parse10x-fnames.1 $inbam` -r "PL:Illumina" $inbam |
    samtools view -bh -o $rgfix0
    samtools reheader <(samtools view -H $rgfix0 | sed '/^@RG\s\+\S\+$/d') $rgfix0 > $rgfix1
    # if order is lexicographic change it to karyotypic
    java -jar $PICARD ReorderSam I=$rgfix1 O=$outbam REFERENCE=$REFSEQ &&
        samtools index $outbam && rm $rgfix0 $rgfix1 $inbam
}

for bam; do
    fix1bam $bam &>> "$destdir/`basename $0`-$starttime.log"
done
