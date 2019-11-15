#! /usr/bin/env python3

import subprocess
import glob
import os
import os.path
import subprocess
import attila_utils

resultdir = '/projects/bsm/attila/results/2019-11-13-panel-of-normals'
NUMBER_THREADS = 16
REFSEQ = '/projects/shared/refgenome/GRCh37/dna/hs37d5.fa'
sentieon_executable = '/opt/sentieon-genomics/current/bin/sentieon'


def bam2pon(bam='/projects/bsm/alignments/PITT_118/PITT_118_NeuN_pl.bam',
        nthreads=NUMBER_THREADS, refseq=REFSEQ, algo='TNhaplotyper'):
    '''
    Create a panel of normal (PON) VCF based on a single BAM

    Arguments
    bam: path to the BAM
    nthreads: number of threads for sentieon driver
    refseq: reference sequence for sentieon driver
    algo: normally TNhaplotyper; hasn't been tested with TNhaplotyper2

    Value:
    a completed process (see subprocess)
    '''
    nthreads = str(nthreads)
    vcf = os.path.basename(bam).replace('.bam', '.vcf')
    outdir = resultdir + os.path.sep + 'VCFs'
    vcf = outdir + os.path.sep + vcf
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    #args = [sentieon_executable, 'driver', '--interval', '22:20000000-25000000', # for testing
    args = [sentieon_executable, 'driver',
            '-t', nthreads, '-r', refseq, '-i', bam, '--algo', algo, '--detect_pon', vcf]
    proc = subprocess.run(args)
    def index_pon_vcf():
        '''
        Replace PON VCF with its gzipped, indexed, version
        '''
        vcfgz = vcf + '.gz'
        args0 = ['bcftools', 'view', '-O', 'z', '-o', vcfgz, vcf]
        proc0 = subprocess.run(args0)
        args1 = ['bcftools', 'index', '--tbi', vcfgz]
        proc1 = subprocess.run(args1)
        if proc1.returncode == 0:
            os.remove(vcf)
            os.remove(vcf + '.idx')
        return(vcfgz)
    val = index_pon_vcf()
    return(val)


def all_bam2pon_merge(bamlist=glob.glob('/projects/bsm/alignments/[MP][SI][ST][MT]_*/*.bam'),
        ponvcf=resultdir + os.path.sep + 'bsm-cmc-pon.vcf.gz'):
    '''
    Create and merge PON VCFs for a list of BAMs

    Arguments:
    bamlist: the list of BAMs (paths)
    ponvcf: path to the single merged PON VCF

    Value:
    a completed process (see subprocess)
    '''
    vcflist = [bam2pon(b) for b in bamlist]
    outdir = resultdir + os.path.sep + 'VCFs'
    vcflist = glob.glob(outdir + os.path.sep + '*.vcf.gz')
    args0 = ['bcftools', 'merge', '-m', 'all', '--force-samples', '-Oz', '-o', ponvcf] + vcflist
    proc0 = subprocess.run(args0)
    args1 = ['bcftools', 'index', '--tbi', ponvcf] + vcflist
    proc1 = subprocess.run(args1, capture_output=True)
    return(proc1)
