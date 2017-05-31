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

for v in snvs indels; do
    for t in {NeuN_mn,muscle}-NeuN_pl; do
        mu2_str_isec $t $v
    done
    reft1vcf=$outmaindir/$subdircaller/NeuN_mn-NeuN_pl/$v/0003.vcf
    reft2vcf=$outmaindir/$subdircaller/muscle-NeuN_pl/$v/0003.vcf
    reftoutdir=$outmaindir/$subdirreftis/$v
    test -d $reftoutdir && rm -r $reftoutdir
    mkdir -p $reftoutdir
    reft1bcf=$reftoutdir/NeuN_mn-NeuN_pl.bcf
    reft2bcf=$reftoutdir/muscle-NeuN_pl.bcf
    vcf2indexedbcf $reft1vcf $reft1bcf $v
    vcf2indexedbcf $reft2vcf $reft2bcf $v
    bcftools isec -p $reftoutdir $reft1bcf $reft2bcf
done

# step 2: summarize results with call set sizes

dosummary () {
    indir=$1
    tmp1=`mktemp`
    tmp2=`mktemp`
    for v in 000{0..2}.vcf; do
        # line numbers with header = set size + 1
        # 'bcftools stats' would be an alternative
        linenowheader=$(grep -v '^##' $indir/$v | wc -l)
        setsize=$(( $linenowheader - 1 ))
        echo $setsize
    done > $tmp1
    readme=$indir/README.txt
    sed -e \
    "1,/^Using/ d;
    s|$indir||g;
    s/for records.*\(private to.*$\|shared by.*$\)/\1/;
    /0003\.vcf/ d" $readme > $tmp2
    paste $tmp1 $tmp2 > $indir/callset-sizes.tsv && rm $tmp1 $tmp2
}

maindir=$HOME/projects/bsm/results/2017-05-29-vcf-comparisons
for t in snvs indels; do
    dosummary $maindir/$filtdir/2_cmp-reftissues/$t/
    for ref in NeuN_mn muscle; do
        dosummary $maindir/$filtdir/1_isec-callers/$ref-NeuN_pl/$t/
    done
done
