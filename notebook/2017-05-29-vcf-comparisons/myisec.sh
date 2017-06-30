#!/usr/bin/env bash

# concordance of call sets under indir from the strelka-mutect2 pilot analysis

usage="usage: ./`basename $0` indir [mutect2_filter]"

# step 1: set operations on call sets

indir=${1:-$HOME/projects/bsm/results/2017-05-03-strelka-mutect2-pilot/32MB}
mutect2_filter=${2}

outmaindir=$HOME/projects/bsm/results/2017-05-29-vcf-comparisons
if test -z $mutect2_filter; then
    filtdir=mutect2-unfilt
else
    filtdir=mutect2-$mutect2_filter
fi
subdircaller=$filtdir/1_isec-callers
subdirreftis=$filtdir/2_cmp-reftissues

# convert VCF into indexed BCF
vcf2indexedbcf () {
    inputf=$1 outputf=$2 vartype=$3 filter=$4
    # note that bcftools uses type 'snps' also for 'snvs'
    if test $vartype == snvs; then
        vartype=snps
    fi
    # filter for vartype, save as .bcf, and index
    if test -z $filter; then
        bcftools view -o $outputf -O b --types $vartype $inputf
    else
        bcftools view -o $outputf -O b --types $vartype -f $filter $inputf
    fi
    bcftools index $outputf
}

# intersection of mutect2 and strelka sets
mu2_str_isec () {
    # args
    tissuepair=$1 # NeuN_mn-NeuN_pl or muscle-NeuN_pl
    vartype=$2 # snvs or indels
    # input
    #indir=$HOME/projects/bsm/results/2017-05-03-strelka-mutect2-pilot/32MB
    inmu2="$indir/$tissuepair-mutect2/out.vcf"
    instr="$indir/$tissuepair-strelka/results/all.somatic.$vartype.vcf"
    # output
    outdir=$outmaindir/$subdircaller/$tissuepair/$vartype
    outmu2=$outdir/mutect2.bcf
    outstr=$outdir/strelka.bcf
    # make outdir
    test -d $outdir && rm -r $outdir
    mkdir -p $outdir
    vcf2indexedbcf $inmu2 $outmu2 $vartype $mutect2_filter
    vcf2indexedbcf $instr $outstr $vartype
    # perform comparison
    bcftools isec -p $outdir $outmu2 $outstr
}

# get number of records in various subsets (intersection, setdifferences)
dosummary1 () {
    inputd=$1
    tmp1=`mktemp`
    tmp2=`mktemp`
    for v in $inputd/000{0..2}.vcf; do
        # line numbers with header = set size + 1
        # 'bcftools stats' would be an alternative
        linenowheader=$(grep -v '^##' $v | wc -l)
        setsize=$(( $linenowheader - 1 ))
        echo $setsize
    done > $tmp1
    readme=$inputd/README.txt
    sed -e \
    "1,/^Using/ d;
    s|$inputd||g;
    s/for records.*\(private to.*$\|shared by.*$\)/\1/;
    /0003\.vcf/ d" $readme > $tmp2 # 0003 has the same no records as 0002 and can be deleted
    paste $tmp1 $tmp2 > $inputd/callset-sizes.tsv && rm $tmp1 $tmp2
}

# outer loop: type of mutation mut
for mut in snvs indels; do
    # clean up subdirectories
    reftoutdir=$outmaindir/$subdirreftis/$mut
    test -d $reftoutdir && rm -r $reftoutdir
    mkdir -p $reftoutdir
    # inner loop: tissue pair t
    for t in NeuN_mn-NeuN_pl muscle-NeuN_pl muscle-NeuN_mn; do
        # pairwise comparison between mutect2 and strelka
        mu2_str_isec $t $mut
        dosummary1 $outmaindir/$subdircaller/$t/$mut
        # prepare intersection of mutect2 and strelka for later comparisons
        tispairvcf=$outmaindir/$subdircaller/$t/$mut/0003.vcf
        tispairbcf=$reftoutdir/$t.bcf
        vcf2indexedbcf $tispairvcf $tispairbcf $mut
    done
    # three-way comparison using the three tissue pairs
    tsv2=$reftoutdir/callset-sizes.tsv # save results in this .tsv file
    tispairlist="muscle-NeuN_pl\tNeuN_mn-NeuN_pl\tmuscle-NeuN_mn"
    echo -e "nrec\tABC\ttissue_pair_A\ttissue_pair_B\ttissue_pair_C" > $tsv2
    # inner loop: subsets defined by bitmap
    # A: muscle-NeuN_pl, B: NeuN_mn-NeuN_pl, C: muscle-NeuN_mn
    # 100: A \ (B ∪ C)
    # 010: B \ (A ∪ C)
    # 001: C \ (A ∪ B)
    # 110: (A ∩ B) \ C
    # 101: (A ∩ C) \ B
    # 011: (B ∩ C) \ A
    # 111: A ∩ B ∩ C
    for bitmap in 100 010 001 110 101 011 111; do
        bcftools isec -n~$bitmap -w 1 -o $reftoutdir/$bitmap.vcf \
            $reftoutdir/{muscle-NeuN_pl.bcf,NeuN_mn-NeuN_pl.bcf,muscle-NeuN_mn.bcf}
        numrecords=$(bcftools stats $reftoutdir/$bitmap.vcf | \
            sed -n '/^.*number of records:\s*\([[:digit:]]\+\)$/ { s//\1/; p }')
        echo -e "$numrecords\t$bitmap\t$tispairlist"
    done >> $tsv2
done
