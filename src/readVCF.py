#! /usr/bin/env python3

import pandas as pd
import numpy as np
import subprocess
import io
import re
import synapseclient
import os.path

def read_annotlist(annotpath='/home/attila/projects/bsm/tables/VCF-HC.annotations', withFORMAT=False):
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
    s = re.sub('^((MSSM|PITT)_([0-9]+))_(.*)$', '\\1:\\4', bsm_sample)
    l = s.split(':')
    indiv_id = 'CMC_' + l[0]
    tissue = l[1]
    return((indiv_id, tissue))

def readVCF(vcfpath, annotlist=read_annotlist()):
    formatstr = make_formatstr(annotlist)
    colnames = [y.replace('INFO/', '') for y in annotlist]
    cmd = ['bcftools', 'query', '-f', formatstr, vcfpath]
    p =  subprocess.run(cmd, capture_output=True)
    df = pd.read_csv(io.BytesIO(p.stdout), sep='\t', names=colnames)
    return(df)

def add_clinical(vcfpath, clinical):
    vcf = readVCF(vcfpath)
    sample = sample_fromVCF(vcfpath)
    indiv_id, tissue = convert_sample(sample)
    df = clinical.loc[[indiv_id]]
    id_df = pd.DataFrame({'Sample': [sample], 'Individual ID': [indiv_id],
        'Tissue': [tissue]})
    id_df.index = df.index
    df = pd.concat([id_df, df], axis=1)
    df = df.iloc[np.zeros(vcf.shape[0]), :]
    df.index = vcf.index
    outvcf = pd.concat([vcf, df], axis=1)
    return(outvcf)

def annotateVCF(invcf='/home/attila/projects/bsm/results/calls/filtered/MSSM_106_brain.ploidy_50.filtered.vcf',
        sample='MSSM_106_NeuN_pl',
        targetdir='/home/attila/projects/bsm/results/calls/annotated/'):
    script = '/home/attila/projects/bsm/src/annotate-vcf-bsm'
    cmd = [script, '-t', targetdir, invcf, sample]
    p =  subprocess.run(cmd, capture_output=True)
    return(p)

def annotate_or_readVCFs(vcflistpath='/big/results/bsm/calls/filtered-vcfs.tsv',
        vcfdir='/home/attila/projects/bsm/results/calls/', fun=annotateVCF, cmc_clinical=None):
    vcflist = pd.read_csv(vcflistpath, sep='\t', names=['sample', 'file'], index_col='sample')
    def helper(sample):
        subdir = 'filtered' if fun is annotateVCF else 'annotated'
        invcf = vcfdir + os.path.sep + subdir + os.path.sep + vcflist.loc[sample, 'file']
        arg2 = sample if cmc_clinical is None else cmc_clinical
        val = fun(invcf, arg2)
        return(val)
    l = [helper(y) for y in vcflist.index]
    return(l)

def annotateVCFs(vcflistpath='/big/results/bsm/calls/filtered-vcfs.tsv',
        vcfdir='/home/attila/projects/bsm/results/calls/'):
    pp = annotate_or_readVCFs(vcflistpath=vcflistpath, vcfdir=vcfdir, fun=annotateVCF, cmc_clinical=None)
    return(pp)

def readVCFs(vcflistpath='/big/results/bsm/calls/filtered-vcfs.tsv',
        vcfdir='/home/attila/projects/bsm/results/calls/'):
    # CMC_Human_clinical_metadata.csv
    syn = synapseclient.login()
    wdir='/tmp'
    cmc_clinical_syn = syn.get('syn2279441', downloadLocation=wdir, ifcollision='overwrite.local')
    cmc_clinical = pd.read_csv(cmc_clinical_syn.path, index_col='Individual ID')
    l = annotate_or_readVCFs(vcflistpath=vcflistpath, vcfdir=vcfdir, fun=add_clinical, cmc_clinical=cmc_clinical)
    val = pd.concat(l, axis=0)
    csvpath = vcfdir + os.path.sep + 'annotations.csv'
    val.to_csv(csvpath, index=False)
    return(val)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', help='main VCF directory (/home/attila/projects/bsm/results/calls/)',
            default='/home/attila/projects/bsm/results/calls/')
    parser.add_argument('-l', '--vcflist', help='list of samples and VCF files (/big/results/bsm/calls/filtered-vcfs.tsv)',
            default='/big/results/bsm/calls/filtered-vcfs.tsv')
    args = parser.parse_args()
    annotateVCFs(vcflistpath=args.vcflist, vcfdir=args.dir)
    readVCFs(vcflistpath=args.vcflist, vcfdir=args.dir)
