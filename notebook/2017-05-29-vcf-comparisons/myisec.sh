#!/usr/bin/env bash

# get intersection from a mutect2-strelka pair of vcf files

outmaindir=$HOME/projects/bsm/results/2017-05-29-vcf-comparisons

vcf2indexedbcf () {
    # filter for vartype, save as .bcf, and index
    bcftools view -o $2 -O b --types $3 $1
    bcftools index $2
}

mu2_str_isec () {
    # args
    tissuepair=$1 # NeuN_mn-NeuN_pl or muscle-NeuN_pl
    vartype=$2 # snps or indels
    # note that the strelka generated vcf file name contains 'snvs' instead of 'snps'
    case $vartype in
        snps) strelkatype=snvs;;
        indels) strelkatype=indels;;
    esac
    # input
    indir=$HOME/projects/bsm/results/2017-05-03-strelka-mutect2-pilot/32MB
    inmu2="$indir/$tissuepair-mutect2/out.vcf"
    instr="$indir/$tissuepair-strelka/results/all.somatic.$strelkatype.vcf"
    # output
    outdir=$outmaindir/caller-isec/$tissuepair/$vartype
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

for v in snps indels; do
    for t in {NeuN_mn,muscle}-NeuN_pl; do
        mu2_str_isec $t $v
    done
    reft1vcf=$outmaindir/caller-isec/NeuN_mn-NeuN_pl/$v/0003.vcf
    reft2vcf=$outmaindir/caller-isec/muscle-NeuN_pl/$v/0003.vcf
    reftoutdir=$outmaindir/cmp-reftissues/$v
    test -d $reftoutdir && rm -r $reftoutdir
    mkdir -p $reftoutdir
    reft1bcf=$reftoutdir/NeuN_mn-NeuN_pl.bcf
    reft2bcf=$reftoutdir/muscle-NeuN_pl.bcf
    vcf2indexedbcf $reft1vcf $reft1bcf $v
    vcf2indexedbcf $reft2vcf $reft2bcf $v
    bcftools isec -p $reftoutdir $reft1bcf $reft2bcf
done
