#! /usr/bin/env python3

import re
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

    Details:
    https://support.sentieon.com/manual/TNseq_usage/tnseq/#generating-a-panel-of-normal-vcf-file
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
    args1 = ['bcftools', 'index', '--tbi', ponvcf]
    proc1 = subprocess.run(args1, capture_output=True)
    return(proc1)


def tnseq_call(bam='/projects/bsm/alignments/PITT_118/PITT_118_NeuN_mn.bam',
        pon='/projects/bsm/attila/results/2019-11-13-panel-of-normals/pon-v1.vcf.gz',
        outdir='/projects/bsm/calls',
        nthreads=NUMBER_THREADS, refseq=REFSEQ, algo='TNhaplotyper'):
    '''
    Call variants with TNseq's TNhaplotyper by default using a PON VCF

    Arguments:
    bam: path to the input BAM
    pon: path to the PON VCF
    outdir: the main directory for the output subdirectories and VCF file
    nthreads: the total number of threads
    refseq: reference genome sequence
    algo: TNhaplotyper by default

    Value:
    a subprocess process object for the sentieon TNseq caller
    '''
    def postproc_vcf(invcf, vcfmaindir, vartype='snps'):
        '''
        Postprocess VCF by filtering for SNPs or indels

        Arguments:
        invcf: input VCF
        vcfmaindir: the main VCF dir
        vartype: snps|indels

        Value:
        a subprocess procedure
        '''
        germ2som = {'snps': 'snvs', 'indels': 'indels'}
        somvt = germ2som[vartype]
        outdir = vcfmaindir + os.sep + somvt
        outvcf = outdir + os.sep + os.path.basename(invcf)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        args = ['bcftools', 'view', '-o', outvcf, '-O', 'z', '-v', vartype, invcf]
        proc = subprocess.run(args, capture_output=True)
        return(proc)
    addthreads = str(nthreads - 1)
    nthreads = str(nthreads)
    # get sample name from BAM header
    args = ['samtools', 'view', '-H', bam]
    proc = subprocess.run(args, capture_output=True)
    headerlines = proc.stdout.decode('utf-8').split('\n')
    matches = [re.match('^@RG.*', line) for line in headerlines]
    samples = set([re.sub('^.*SM:([^\t]+)\t.*$', '\\1', m.string) for m in
        matches if m])
    if len(samples) != 1:
        raise Exception('There must be one and only one sample in the BAM header.  Quitting...')
    sample = samples.pop()
    # get subject and tissue
    # this works with MSSM373_NeuN_pl, MSSM_373_NeuN_pl or Common_7_NeuN_pl
    pattern = '([a-zA-Z]+)_?([0-9]+)_(.*)$'
    subject = re.sub(pattern, '\\1_\\2', sample)
    tissue = re.sub(pattern, '\\3', sample)
    bname = subject + '_' + tissue
    if bname != os.path.basename(bam).replace('.bam', ''):
        raise Exception('Sample name from BAM header does not match that from filename.  Quitting...')
    # output dir and VCF
    vcfmaindir = outdir + os.path.sep + subject + os.path.sep + tissue  + '-PON'
    vcfdir = vcfmaindir + os.path.sep + algo
    if not os.path.exists(vcfdir):
        os.makedirs(vcfdir)
    vcf = vcfdir + os.path.sep + algo + '.vcf.gz'
    newpon = pon_without_sample(bam=bam, addthreads=addthreads, mergedpon=pon,
            subtract=False) # remove sample from PON
    #args = [sentieon_executable, 'driver', '--interval', '22:20000000-25000000', # for testing
    args = [sentieon_executable, 'driver',
            '-t', nthreads, '-r', refseq, '-i', bam, '--algo', algo,
            '--tumor_sample', sample, '--pon', newpon, vcf]
    proc = subprocess.run(args, capture_output=True)
    # STDERR to logfile
    log = vcf.replace('.vcf.gz', '.log')
    print('Callset:\n' + vcf)
    print('Logfile:\n' + log)
    with open(log, mode='w') as f:
        print(proc.stderr.decode('utf-8'), file=f)
    p_snps = postproc_vcf(invcf=vcf, vcfmaindir=vcfmaindir, vartype='snps')
    p_indels = postproc_vcf(invcf=vcf, vcfmaindir=vcfmaindir, vartype='indels')
    return(proc)


def pon_without_sample(bam, addthreads,
        mergedpon='/projects/bsm/attila/results/2019-11-13-panel-of-normals/pon-v1.vcf.gz',
        subtract=False):
    '''
    Remove sample specific records from PON VCF if necessary

    Arguments:
    bam: path to bam corresponding to the sample
    addthreads: additional threads
    mergedpon: path to the PON VCF to remove the sample from
    subtract: set subtraction: mergedpon - samplepon

    Value:
    path to the PON VCF without the sample

    Details:
    If PON VCF "mergedpon" does not contain sample then no new PON VCF is created;
    in this case the path to the original PON VCF is returned.
    If subtract is True then the PON corresponding to the input bam file is
    subtracted from the merged PON: mergedpon - samplepon.  Otherwise, by
    default, a new merged PON is made excluding the sample PON but including
    all other individual PONs in the
    '/projects/bsm/attila/results/2019-11-13-panel-of-normals/VCFs/'
    directory.
    '''
    addthreads = str(addthreads)
    pondir = os.path.dirname(mergedpon)
    ponbn = os.path.basename(mergedpon).replace('.vcf.gz', '')
    sample = os.path.basename(bam).replace('.bam', '')
    samplepon = pondir + os.sep + 'VCFs' + os.sep + sample + '.vcf.gz'
    newpon = pondir + os.path.sep + ponbn + '-' + sample + '.vcf.gz'
    if not os.path.exists(samplepon):
        return(mergedpon)
    if subtract:
        args = ['bcftools', 'isec', '--complement', '-Oz', '-o', newpon, '-w1',
                '--threads', addthreads, mergedpon, samplepon]
        proc = subprocess.run(args, capture_output=True)
        if proc.returncode == 0:
            args1 = ['bcftools', 'index', '--tbi', newpon]
            proc1 = subprocess.run(args1)
            return(newpon)
        else:
            raise Exception('Unidentified bcftools error.  Quitting...')
    else:
        vcfdir = pondir + os.sep + 'VCFs'
        vcflist = glob.glob(vcfdir + os.sep + '*.vcf.gz')
        # filter for canonical filenames
        vcflist = [y for y in vcflist if re.match(vcfdir + os.sep +
            '(MSSM|PITT)_[0-9]{3}_(muscle|NeuN_pl|NeuN_mn)\.vcf.gz', y)]
        # remove sample PON from the list
        vcflist.remove(samplepon)
        args = ['bcftools', 'merge', '-m', 'all', '--force-samples', '-Oz',
                '-o', newpon, '--threads', addthreads] + vcflist
        proc = subprocess.run(args, capture_output=True)
        return(newpon)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='input BAM file')
    parser.add_argument('-p', '--pon', help='PON (panel of normal) VCF',
            default='/projects/bsm/attila/results/2019-11-13-panel-of-normals/pon-v1.vcf.gz')
    parser.add_argument('-o', '--outdir', help='main output directory',
            default='/projects/bsm/calls')
    parser.add_argument('-t', '--nthreads', help='number of threads',
            default=16, type=int)
    parser.add_argument('-r', '--refseq', help='reference genome sequence',
            default='/projects/shared/refgenome/GRCh37/dna/hs37d5.fa')
    parser.add_argument('-a', '--algo', help='TNseq calling algorithm',
            default='TNhaplotyper')
    args = parser.parse_args()
    proc = tnseq_call(bam=args.bam, pon=args.pon, outdir=args.outdir, nthreads=args.nthreads,
            refseq=args.refseq, algo=args.algo)
    print(proc.stdout.decode('utf-8'))
