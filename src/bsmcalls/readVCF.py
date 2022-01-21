#! /usr/bin/env python3

import pandas as pd
import numpy as np
import subprocess
import io
import re
import os.path
import bsmutils

# Source:
# https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/browserlabelmap_15_coreMarks.tab
state15label = {1 :  "TssA",
                2 :  "TssAFlnk",
                3 :  "TxFlnk",
                4 :  "Tx",
                5 :  "TxWk",
                6 :  "EnhG",
                7 :  "Enh",
                8 :  "ZNF/Rpts",
                9 :  "Het",
                10:  "TssBiv",
                11:  "BivFlnk",
                12:  "EnhBiv",
                13:  "ReprPC",
                14:  "ReprPCWk",
                15:  "Quies"}

'''
========================================================================
Reading VCFs
========================================================================
'''

def read_annotlist(annotpath=bsmutils.get_bsmdir() + '/tables/VCF-HC.annotations', withFORMAT=False):
    '''
    Reads a file containing list of annotations in VCFs into a list.

    Parameters
    annotpath: the path to the aforementioned file
    withFORMAT: if False (default) the FORMAT fields are omitted

    Value: the list of annotations
    '''
    with open(annotpath) as f:
        l = f.readlines()
    l = [y.replace('\n', '') for y in l] # remove newline characters
    if not withFORMAT:
        l = [y for y in l if not re.match('^FORMAT', y)]
    return(l)

def make_formatstr(annotlist):
    formatstr = '%' + '\t%'.join(annotlist) + '\n'
    return(formatstr)

def sample_fromVCF(vcfpath):
    cmd = ['bcftools', 'view', '-h', vcfpath]
    p =  subprocess.run(cmd, capture_output=True)
    colnames = p.stdout.splitlines()[-1]
    sample = re.sub('^.*\t([^\t]+)$', '\\1', colnames.decode('utf-8'))
    return(sample)

def convert_sample(bsm_sample):
    cmc_pattern = '^((MSSM|PITT)_([0-9]+))_(.*)$'
    walsh_pattern = '^([A-Z]+[^_]*)$'
    if re.match(cmc_pattern, bsm_sample):
        s = re.sub(cmc_pattern, '\\1:\\4', bsm_sample)
        l = s.split(':')
        indiv_id = 'CMC_' + l[0]
        tissue = l[1]
    elif re.match(walsh_pattern, bsm_sample):
        indiv_id = bsm_sample
        tissue = 'frontal cortex'
    else:
        raise ValueError('Unexpected sample name: ' + bsm_sample)
    return((indiv_id, tissue))

def readVCF(vcfpath, annotlist=read_annotlist()):
    '''
    Reads the calls/records of VCF into rows of a DataFrame

    Arguments
    vcfpath: path to VCF file
    annotlist: the list of annotations (CHROM, POS, ..., FILTER/1, ..., INFO/1, ...)

    Value:
    calls: a pandas DataFrame
    '''
    formatstr = make_formatstr(annotlist)
    colnames = [y.replace('INFO/', '') for y in annotlist]
    cmd = ['bcftools', 'query', '-f', formatstr, vcfpath]
    p =  subprocess.run(cmd, capture_output=True)
    calls = pd.read_csv(io.BytesIO(p.stdout), sep='\t', names=colnames,
                        na_values='.', dtype={'CHROM': 'object'})
    # extra columns
    l = [state15label[x] for x in calls['ChromatinState_DLPFC']]
    calls['ChromatinState_DLPFC'] = pd.Categorical(l, categories=state15label.values(), ordered=True)
    # SiPhy turned out to refer to hg18/GRCh36 instead of hg19/GRCh37
    #calls['evolConstrain'] = [not np.isnan(y) for y in calls['SiPhyLOD']]
    sample = sample_fromVCF(vcfpath)
    indiv_id, tissue = convert_sample(sample)
    calls['CHROM'] = calls['CHROM']
    calls['Individual ID'] = indiv_id
    calls['Tissue'] = tissue
    calls['Mutation'] = [str(a) + '/' + str(b) for a, b in zip(calls['REF'], calls['ALT'])]
    # set index
    calls = calls.set_index(['Individual ID', 'Tissue', 'CHROM', 'POS', 'Mutation'], drop=True)
    #calls.columns = pd.MultiIndex.from_product([['VCF'], calls.columns], names=['Source', 'Annotation'])
    return(calls)

def readVCFs(vcflistpath=bsmutils.get_bsmdir() + '/results/calls/filtered-vcfs.tsv',
        vcfdir=bsmutils.get_bsmdir() + '/results/calls/', clean=True, annotlist=read_annotlist()):
    '''
    Reads the calls/records of several VCFs into rows of a single DataFrame

    Arguments
    vcflistpath: path to file listing all VCFs
    vcfdir: the directory of the VCFs
    clean: weather to remove redundant & degenerate columns

    Value:
    calls: a pandas DataFrame
    '''
    vcflist = pd.read_csv(vcflistpath, sep='\t', names=['sample', 'file'], index_col='sample')
    vcflist['filepath'] = [vcfdir + os.sep + 'annotated' + os.sep + f for f in vcflist['file']]
    l = [readVCF(y, annotlist=annotlist) for y in vcflist['filepath']]
    calls = pd.concat(l, axis=0)
    if clean:
        calls = clean_calls(calls, dropna=True, dropdegenerate=True, dropredundant=True)
    return(calls)


'''
========================================================================
Annotate VCFs
========================================================================
'''

def clean_calls(calls, dropna=True, dropdegenerate=True, dropredundant=True):
    '''
    Remove redundant & degenerate columns and those with missing values

    Details:
    The set of redundant and degenerate variables was established in 2020-08-07-cleaning-data
    '''
    redundant_vars = ['FILTER/HC', 'FILTER/EXT', 'QUAL', 'AC', 'MLEAC', 'MLEAF', 'QD']
    degenerate_vars = ['AN', 'MQ', 'MQRankSum']
    if dropna:
        calls = calls.dropna(axis=1)
    if dropdegenerate:
        calls = calls.drop(columns=degenerate_vars)
    if dropredundant:
        calls = calls.drop(columns=redundant_vars)
    return(calls)


def annotateVCF(invcf=bsmutils.get_bsmdir() + '/results/calls/filtered/MSSM_106_brain.ploidy_50.filtered.vcf',
        sample='MSSM_106_NeuN_pl',
        targetdir=bsmutils.get_bsmdir() + '/results/calls/annotated/'):
    '''
    Some help would be nice
    '''
    script = bsmutils.get_bsmdir() + '/src/annotate-vcf-bsm'
    cmd = [script, '-t', targetdir, invcf, sample]
    p =  subprocess.run(cmd, capture_output=True)
    return(p)

def annotateVCFs(vcflistpath=bsmutils.get_bsmdir() + '/results/calls/filtered-vcfs.tsv',
        vcfdir=bsmutils.get_bsmdir() + '/results/calls/'):
    vcflist = pd.read_csv(vcflistpath, sep='\t', names=['sample', 'file'], index_col='sample')
    def helper(sample):
        invcf = vcfdir + os.path.sep + 'filtered' + os.path.sep + vcflist.loc[sample, 'file']
        targetdir = vcfdir + os.path.sep + 'annotated' + os.path.sep
        val = annotateVCF(invcf=invcf, sample=sample, targetdir=targetdir)
        return(val)
    pp = [helper(y) for y in vcflist.index]
    return(pp)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', help='main VCF directory (bsm/results/calls/)',
            default=bsmutils.get_bsmdir() + '/results/calls/')
    parser.add_argument('-l', '--vcflist', help='list of samples and VCF files (bsm/results/calls/filtered-vcfs.tsv)',
            default=bsmutils.get_bsmdir() + '/results/calls/filtered-vcfs.tsv')
    args = parser.parse_args()
    annotateVCFs(vcflistpath=args.vcflist, vcfdir=args.dir)
    readVCFs(vcflistpath=args.vcflist, vcfdir=args.dir)
