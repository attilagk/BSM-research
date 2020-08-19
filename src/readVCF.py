#! /usr/bin/env python3

import pandas as pd
import numpy as np
import subprocess
import io
import re
import synapseclient
import os.path

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

def read_annotlist(annotpath='/home/attila/projects/bsm/tables/VCF-HC.annotations', withFORMAT=False):
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
    calls = pd.read_csv(io.BytesIO(p.stdout), sep='\t', names=colnames, na_values='.')
    # extra columns
    calls['ChromatinState_DLPFC'] = [state15label[x] for x in calls['ChromatinState_DLPFC']]
    calls['evolConstrain'] = [not np.isnan(y) for y in calls['SiPhyLOD']]
    sample = sample_fromVCF(vcfpath)
    indiv_id, tissue = convert_sample(sample)
    calls['Individual ID'] = indiv_id
    calls['Tissue'] = tissue
    calls['Mutation'] = [str(a) + '->' + str(b) for a, b in zip(calls['REF'], calls['ALT'])]
    # set index
    calls = calls.set_index(['Individual ID', 'Tissue', 'CHROM', 'POS', 'Mutation'], drop=False)
    return(calls)


def annotateVCF(invcf='/home/attila/projects/bsm/results/calls/filtered/MSSM_106_brain.ploidy_50.filtered.vcf',
        sample='MSSM_106_NeuN_pl',
        targetdir='/home/attila/projects/bsm/results/calls/annotated/'):
    '''
    Some help would be nice
    '''
    script = '/home/attila/projects/bsm/src/annotate-vcf-bsm'
    cmd = [script, '-t', targetdir, invcf, sample]
    p =  subprocess.run(cmd, capture_output=True)
    return(p)

def annotateVCFs(vcflistpath='/big/results/bsm/calls/filtered-vcfs.tsv',
        vcfdir='/home/attila/projects/bsm/results/calls/'):
    vcflist = pd.read_csv(vcflistpath, sep='\t', names=['sample', 'file'], index_col='sample')
    def helper(sample):
        invcf = vcfdir + os.path.sep + 'filtered' + os.path.sep + vcflist.loc[sample, 'file']
        val = annotateVCF(invcf=invcf, sample=sample)
        return(val)
    pp = [helper(y) for y in vcflist.index]
    return(pp)

def readVCFs(vcflistpath='/big/results/bsm/calls/filtered-vcfs.tsv',
        vcfdir='/home/attila/projects/bsm/results/calls/'):
    # CMC_Human_clinical_metadata.csv
    syn = synapseclient.login()
    wdir = '/home/attila/projects/bsm/resources/'
    clinical_syn = syn.get('syn2279441', downloadLocation=wdir, ifcollision='overwrite.local')
    clinical = pd.read_csv(clinical_syn.path, index_col='Individual ID')
    vcflist = pd.read_csv(vcflistpath, sep='\t', names=['sample', 'file'], index_col='sample')
    vcflist['filepath'] = [vcfdir + '/annotated/' + f for f in vcflist['file']]
    l = [readVCF(y) for y in vcflist['filepath']]
    val = pd.concat(l, axis=0)
    return(val)
#    val = set_dtypes(val)
#    val = add_ancestry(val)
#    csvpath = vcfdir + os.path.sep + 'annotations.csv'
#    val.to_csv(csvpath, index=False)
#    val = index_with_variants(val)
#    return(val)

#def set_dtypes(vcf):
#    vcf['Dx'] = pd.Categorical(vcf['Dx'], categories=['Control', 'SCZ'])
#    vcf['chromatinState_DLPFC'] = pd.Categorical(vcf['chromatinState_DLPFC'], categories=state15label.values())
#    return(vcf)

def add_ancestry(callsdf,
        ancestrypath='/home/attila/projects/bsm/resources/cmc-ancestry/CMC_MSSM-Penn-Pitt_DNA_GENOTYPE_ANCESTRY_GemTools.tsv'):
    ancestry = pd.read_csv(ancestrypath, sep='\t', index_col='Individual_ID')
    # take care of missing values
    c = set(callsdf['Individual ID'])
    a = set(ancestry.index)
    ancestry = ancestry.reindex(list(ancestry.index) + list(c - a))
    # reshape ancestry according to callsdf
    ancestry_reshaped = ancestry.loc[callsdf['Individual ID'], :]
    ancestry_reshaped.index = callsdf.index
    val = pd.concat([callsdf, ancestry_reshaped], axis=1)
    return(val)

def agg2samples(calls):
    callsg = calls.groupby('Sample')
    # call specific numeric variables are aggregated by mean and std
    A = callsg.aggregate(
            QUAL_mean=pd.NamedAgg(column='QUAL', aggfunc=np.mean),
            QUAL_std=pd.NamedAgg(column='QUAL', aggfunc=np.std),
            AF_mean=pd.NamedAgg(column='AF', aggfunc=np.mean),
            AF_std=pd.NamedAgg(column='AF', aggfunc=np.std),
            BaseQRankSum_mean=pd.NamedAgg(column='BaseQRankSum', aggfunc=np.mean),
            BaseQRankSum_std=pd.NamedAgg(column='BaseQRankSum', aggfunc=np.std),
            DP_mean=pd.NamedAgg(column='DP', aggfunc=np.mean),
            DP_std=pd.NamedAgg(column='DP', aggfunc=np.std),
            FS_mean=pd.NamedAgg(column='FS', aggfunc=np.mean),
            FS_std=pd.NamedAgg(column='FS', aggfunc=np.std),
            QD_mean=pd.NamedAgg(column='QD', aggfunc=np.mean),
            QD_std=pd.NamedAgg(column='QD', aggfunc=np.std),
            ReadPosRankSum_mean=pd.NamedAgg(column='ReadPosRankSum', aggfunc=np.mean),
            ReadPosRankSum_std=pd.NamedAgg(column='ReadPosRankSum', aggfunc=np.std),
            SOR_mean=pd.NamedAgg(column='SOR', aggfunc=np.mean),
            SOR_std=pd.NamedAgg(column='SOR', aggfunc=np.std),
            VQSLOD_mean=pd.NamedAgg(column='VQSLOD', aggfunc=np.mean),
            VQSLOD_std=pd.NamedAgg(column='VQSLOD', aggfunc=np.std),
            # call specific categorical variables are aggregated by count
            evolConstrain_count=pd.NamedAgg(column='evolConstrain', aggfunc=np.sum),
            # TODO: chromatinState_DLPFC value counts to separate columns for each value like Quies
            )
    # sample specific variables are left as it is
    B = calls.loc[:, 'Sample':].groupby('Sample').aggregate('first')
    samples = pd.concat([A, B], axis=1)
    return(samples)

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
