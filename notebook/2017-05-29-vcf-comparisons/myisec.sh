#!/usr/bin/env bash

# concordance of call sets under indir from the strelka-mutect2 pilot analysis

usage="usage: ./`basename $0` indir"

# step 1: set operations on call sets

indir=${1:-$HOME/projects/bsm/results/2017-05-03-strelka-mutect2-pilot/32MB}
outmaindir=$HOME/projects/bsm/results/2017-05-29-vcf-comparisons
subdircaller=1_isec-callers
subdirreftis=2_cmp-reftissues

vcf2indexedbcf () {
    # note that bcftools uses type 'snps' also for 'snvs'
    vartype=$3
    if test $vartype == snvs; then
        vartype=snps
    fi
    # filter for vartype, save as .bcf, and index
    bcftools view -o $2 -O b --types $vartype $1
    bcftools index $2
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
    vcf2indexedbcf $inmu2 $outmu2 $vartype
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
    dosummary $maindir/2_cmp-reftissues/$t/
    for ref in NeuN_mn muscle; do
        dosummary $maindir/1_isec-callers/$ref-NeuN_pl/$t/
    done
done
