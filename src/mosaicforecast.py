#! /usr/bin/env python3

import subprocess
import os
import os.path
import shutil
import pandas as pd
import re

SegDup_and_clustered_bed = '/home/attila/projects/MosaicForecast/resources/SegDup_and_clustered.bed'
gnomAD_genome_VCF = '/projects/shared/gnomAD/gnomad.exomes.r2.1.1.sites.1.vcf.bgz'
yifans_filter_script = '/home/attila/projects/MosaicForecast/MuTect2-PoN_filter.py'
NUMBER_THREADS=16


def MF_recommended_filter_yifan(invcf, nthreads=NUMBER_THREADS, keepVCF=False, mt_pon_filter_name='filt_segdup_clust'):
    '''
    Wrapper around Yifan's MuTect2-PoN_filter.py

    Arguments
    invcf: path to input VCF
    nthreads: number of threads
    keepVCF: should we keep the temporary VCF? (only matters if invcf is gzipped)
    mt_pon_filter_name: the name of the subdirectory where the output VCF is written to

    Value
    a subprocess process object representing the last filtering step

    Details
    Yifan's filter script is extremely inefficient both in terms of time and
    memory.  The memory requirement is so large that for a typical VCF for our
    BSM study not only 32 GiB physical memory but also the swap memory is
    filled up.  The inefficiency follows from the fact that the script
    performs the filtering steps by reading all data in memory at once and
    then uses Python's re module to find VCF records that match the desired
    filters.  This made me reimplement Yifan's filter script using the efficient
    bcftools suite.
    '''
    def filter_vcf_for_bed(invcf, outvcf, bed, nthreads=NUMBER_THREADS):
        addthreads = str(nthreads - 1)
        regionsf = bed2regions_file(bed)
        args = ['bcftools', 'view', '--threads', addthreads, '-R', regionsf, '-Oz', '-o', outvcf, invcf]
        proc = subprocess.run(args, capture_output=True)
        return(proc)
    addthreads = str(nthreads - 1)
    # infer filetype
    gzmatch = re.match('.*\.vcf.gz$', invcf)
    if gzmatch:
        invcfgz = invcf
        vcfext = '.vcf.gz'
        tmpvcf = invcf.replace('.gz', '')
        args = ['bcftools', 'view', '--threads', addthreads, '-Ov', '-o', tmpvcf, invcf]
        proc = subprocess.run(args, capture_output=True)
        pass
    elif re.match('.*\.vcf$', invcf):
        invcfgz = invcf + '.gz'
        args = ['bcftools', 'view', '--threads', addthreads, '-Oz', '-o', invcfgz, invcf]
        proc = subprocess.run(args, capture_output=True)
        args = ['bcftools', 'index', '--tbi', '--threads', addthreads, invcfgz]
        proc = subprocess.run(args, capture_output=True)
        vcfext = '.vcf'
        tmpvcf = invcf
    else:
        raise Exception('Error: ' + invcf + ' is not a .vcf or .vcf.gz file')
    # directory and file names
    vcfbname = os.path.basename(invcf).replace(vcfext, '')
    tmpbed = vcfbname + '.bed'
    filter_dir = os.path.dirname(invcf) + os.sep + mt_pon_filter_name
    outbed = filter_dir + os.sep + vcfbname + '.bed'
    outvcf = filter_dir + os.sep + vcfbname + '.vcf.gz'
    os.makedirs(filter_dir, exist_ok=True)
    # get bed with positions to keep
    args = ['python3', yifans_filter_script, vcfbname, tmpvcf, SegDup_and_clustered_bed]
    proc = subprocess.run(args, capture_output=True)
    shutil.move(tmpbed, outbed)
    # perform filtering
    proc = filter_vcf_for_bed(invcf=invcfgz, outvcf=outvcf, bed=outbed,
            nthreads=nthreads)
    # clean up
    if gzmatch and not keepVCF:
        os.remove(tmpvcf)
    return(proc)


def bed2regions_file(bedfile):
    # this depends on the MuTect2-PoN_filter.py script
    colnames = ['chr', 'pos0', 'pos1', 'ref', 'alt', 'sample', 'depth', 'AF']
    bed = pd.read_csv(bedfile, sep='\t', header=None, names=colnames)
    # positions in bed files are 0 based so add 1 to get 1 based positions
    bed['pos'] = bed['pos0'] + 1
    regfname = os.path.dirname(bedfile) + os.sep + os.path.basename(bedfile).replace('.bed', '.regions')
    reg = bed.loc[:, ['chr', 'pos']]
    reg.to_csv(regfname, sep='\t', header=False, index=False)
    return(regfname)


def prefilter(invcf, outvcf):
    '''
    Perform first filtering steps to remove certain calls (see details).

    Arguments
    invcf: input VCF as .vcf.gz
    outvcf: output VCF as .vcf.gz

    Value:
    the process object corresponding to the filtering step

    Details
    The following calls are removed:
    FILTER: panel_of_normals, str_contraction, triallelic_site, t_lod_fstar
    AD[0:1] < 2 (depth of alt allele)
    AF > 0.4 (allele frequency)
    AF < 0.02 if FORMAT/PGT is 0|1, where PGT is Pyhysical phasing haplotype info
    AF < 0.02 otherwise (originally this was AF < 0.03 but changed to match behavior of Yifan's filter script)
    '''
    filtvalues = ['"panel_of_normals"', '"str_contraction"', '"triallelic_site"', '"t_lod_fstar"'] 
    expr1_l = ['FILTER!=' + s for s in filtvalues]
    expr2_l = ['AF<=0.4', 'AD[0:1]>=2' ]
    expr3A_l = ['AF>=0.02', 'FORMAT/PGT==0|1']
    #expr3B_l = ['AF>0.03', 'FORMAT/PGT!=0|1']
    expr3B_l = ['AF>=0.02', 'FORMAT/PGT!=0|1']
    exprA = ' && '.join(expr1_l + expr2_l + expr3A_l)
    exprB = ' && '.join(expr1_l + expr2_l + expr3B_l)
    expr = exprA + ' || ' + exprB
    args = ['bcftools', 'view', '--include', expr, '-Oz', '-o', outvcf, invcf]
    proc = subprocess.run(args, capture_output=True)
    args1 = ['bcftools', 'index', '--tbi', outvcf]
    proc1 = subprocess.run(args1, capture_output=True)
    return(proc)


def segdup_clustered_filter(invcf, replaceinvcf=False):
    '''
    Remove calls in segmental duplications and clutered regions

    Arguments
    invcf: pathname to input vcf.gz file
    replaceinvcf: whether to replace invcf with outvcf

    Value:
    the last process object of subprocess module

    Details:
    If replaceinvcf is True then the output .vcf.gz file replaces invcf but
    only after a copy of invcf is saved with the -pre.vcf.gz suffix.
    Otherwise the name of invcf stays the same and a
    '-segdup_clustered.vcf.gz' suffix is added to the output VCF.
    '''
    if replaceinvcf:
        outvcf = invcf
        # save the original invcf (and its index) with the -pre.vcf.gz suffix
        invcf = invcf.replace('.vcf.gz', '-pre.vcf.gz')
        shutil.move(outvcf, invcf)
        shutil.move(outvcf + '.tbi', invcf + '.tbi')
    else:
        # the ouput VCF has the '-segdup_clustered.vcf.gz' suffix
        outvcf = os.path.dirname(invcf) + os.sep + os.path.basename(invcf).replace('.vcf.gz', '') + '-segdup_clustered' + '.vcf.gz'
    inbed = invcf.replace('.vcf.gz', '.bed')
    outbed = inbed.replace('.bed', '-segdup_clustered.bed')
    # we need non-gzipped VCF for vcf2bed conversion
    args = ['bcftools', 'view', '-Ov', invcf]
    proc = subprocess.run(args, capture_output=True)
    # using bedops convert VCF to BED for bedtools arithmetic
    # https://bedops.readthedocs.io/en/latest/index.html
    with open(inbed, 'w') as f:
        proc1 = subprocess.run(['vcf2bed'], input=proc.stdout, stdout=f)
    # get regions in BED format to retain calls in
    # this uses subtractBed a.k.a bedtools subtract
    args2 = ['subtractBed', '-a', inbed, '-b', SegDup_and_clustered_bed]
    with open(outbed, 'w') as f:
        proc2 = subprocess.run(args2, stdout=f)
    # retain calls in the regions
    args3 = ['bcftools', 'view', '-Oz', '-o', outvcf, '--regions-file', outbed, invcf]
    proc3 = subprocess.run(args3, capture_output=True)
    # index output VCF
    args4 = ['bcftools', 'index', '--tbi', outvcf]
    proc4 = subprocess.run(args4, capture_output=True)
    # cleanup
    os.remove(inbed)
    os.remove(outbed)
    if replaceinvcf:
        os.remove(invcf)
        os.remove(invcf + '.tbi')
    return(proc3)


def bcftools_pipe(cmd, invcf, outvcf=None):
    '''
    Extends and runs bcftools cmd depending on if invcf and/or outvcf is None

    Arguments
    cmd: the argument list of the bcftools cmd e.g ['bcftools', 'annotate',...]
    invcf: a subprocess.CompletedProcess object if STDIN or a string evaluating to the path to the input VCF
    outvcf: None if STDOUT or a string evaluating to the path to the output VCF

    Value
    a subprocess.CompletedProcess object
    '''
    # Is output STDOUT?
    if outvcf is not None:
        cmd += ['-o', outvcf]
    # Is input STDIN?
    if type(invcf) is subprocess.CompletedProcess:
        args = cmd + ['-']
        proc = subprocess.run(args, input=invcf.stdout, capture_output=True)
    else:
        args = cmd + [invcf]
        proc = subprocess.run(args, capture_output=True)
    return(proc)


def gnomAD_AF_annotate(invcf, outvcf, nthreads=NUMBER_THREADS, gnomADvcf=gnomAD_genome_VCF):
    # annotate with gnomAD allele frequency
    addthreads = str(nthreads - 1)
    cmd = ['bcftools', 'annotate', '--threads', addthreads, '-a', gnomADvcf, '-c', 'INFO/AF', '-Oz']
    proc = bcftools_pipe(cmd=cmd, invcf=invcf, outvcf=outvcf)
    #proc = subprocess.run(args, capture_output=True)
    return(proc)

def gnomAD_AF_filter(invcf, outvcf, AFthrs=0.001):
    AFthrs = str(AFthrs)
    args = ['bcftools', 'view', '--exclude', 'INFO/AF>=' + AFthrs, '-Oz', '-o', outvcf, invcf]
    proc = subprocess.run(args, capture_output=True)
    return(proc)


def MF_recommended_filter(invcf, do_gnomAD=True, subdir='MosaicForecast_input'):
    '''
    Perform all filtering steps on the MuTect2-PON callset
    TODO: gnomAD filter
    '''
    filter_dir = os.path.dirname(invcf) + os.sep + 'filtered' + os.sep + subdir
    os.makedirs(filter_dir, exist_ok=True)
    outvcf = filter_dir + os.sep + os.path.basename(invcf)
    prefilter(invcf=invcf, outvcf=outvcf)
    proc = segdup_clustered_filter(invcf=outvcf, replaceinvcf=True)
    return(proc)
