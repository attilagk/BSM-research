#! /usr/bin/env python3

import subprocess
import os
import os.path
import shutil
import pandas as pd
import re

SegDup_and_clustered_bed = '/home/attila/projects/MosaicForecast/resources/SegDup_and_clustered.bed'
mt_pon_filter_script = '/home/attila/projects/MosaicForecast/MuTect2-PoN_filter.py'
mt_pon_filter_name = 'filt_segdup_clust'
NUMBER_THREADS=16

def mt_pon_filter(invcf, nthreads=NUMBER_THREADS, keepVCF=False):
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
    args = ['python3', mt_pon_filter_script, vcfbname, tmpvcf, SegDup_and_clustered_bed]
    proc = subprocess.run(args, capture_output=True)
    shutil.move(tmpbed, outbed)
    # perform filtering
    proc = filter_vcf_for_bed(invcf=invcfgz, outvcf=outvcf, bed=outbed,
            nthreads=nthreads)
    # clean up
    if gzmatch and not keepVCF:
        os.remove(tmpvcf)
    return(proc)


def my_prefilter(invcf, outvcf):
    filtvalues = ['"panel_of_normals"', '"str_contraction"', '"triallelic_site"', '"t_lod_fstar"'] 
    expr1_l = ['FILTER!=' + s for s in filtvalues]
    expr2_l = ['AF<0.4', 'AD[0:1]>2' ]
    expr3A_l = ['AF>0.02', 'FORMAT/PGT==0|1']
    expr3B_l = ['AF>0.03', 'FORMAT/PGT!=0|1']
    exprA = ' && '.join(expr1_l + expr2_l + expr3A_l)
    exprB = ' && '.join(expr1_l + expr2_l + expr3B_l)
    expr = exprA + ' || ' + exprB
    args = ['bcftools', 'view', '--include', expr, '-Oz', '-o', outvcf, invcf]
    proc = subprocess.run(args, capture_output=True)
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
