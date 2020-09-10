#! /usr/bin/env python3

import pandas as pd
import numpy as np
import subprocess
import io
import re
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

'''
========================================================================
Reading VCFs
========================================================================
'''

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
    calls = pd.read_csv(io.BytesIO(p.stdout), sep='\t', names=colnames, na_values='.')
    # extra columns
    l = [state15label[x] for x in calls['ChromatinState_DLPFC']]
    calls['ChromatinState_DLPFC'] = pd.Categorical(l, categories=state15label.values(), ordered=True)
    calls['evolConstrain'] = [not np.isnan(y) for y in calls['SiPhyLOD']]
    sample = sample_fromVCF(vcfpath)
    indiv_id, tissue = convert_sample(sample)
    calls['Individual ID'] = indiv_id
    calls['Tissue'] = tissue
    calls['Mutation'] = [str(a) + '/' + str(b) for a, b in zip(calls['REF'], calls['ALT'])]
    # set index
    calls = calls.set_index(['Individual ID', 'Tissue', 'CHROM', 'POS', 'Mutation'], drop=True)
    return(calls)

def readVCFs(vcflistpath='/big/results/bsm/calls/filtered-vcfs.tsv',
        vcfdir='/home/attila/projects/bsm/results/calls/annotated/', clean=True):
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
    vcflist['filepath'] = [vcfdir + f for f in vcflist['file']]
    l = [readVCF(y) for y in vcflist['filepath']]
    calls = pd.concat(l, axis=0)
    if clean:
        calls = clean_calls(calls, dropna=True, dropdegenerate=True, dropredundant=True)
    return(calls)

'''
========================================================================
SNPnexus annotations
========================================================================
'''

def read_TXT_per_annotation(tsvpath, indivID, tissue='NeuN_pl'):
    '''
    Reads a TXT_per_annotation file of SNPnexus into an indexed DataFrame
    '''
    annot = pd.read_csv(tsvpath, sep='\t')
    def varid2index(varid):
        s = re.sub('^chr(.+):1$', '\\1', varid)
        val = s.split(':')
        return(val)
    l = [[indivID, tissue] + varid2index(s) for s in annot['Variation ID']]
    a = np.array(l)
    columns = ['Individual ID', 'Tissue', 'CHROM', 'POS', 'Mutation']
    df = pd.DataFrame(a, columns=columns)
    df['POS'] = df['POS'].astype('int64')
    annot.index = pd.MultiIndex.from_frame(df)
    return(annot)

def annotation_duplicates(annot, sep=':'):
    '''
    Takes care of rows with duplicated index that occur e.g with overlapping genes

    Arguments
    annot: a pandas DataFrame with possible duplicates
    sep: the separator for the collapsed list of strings

    Value: a copy of annot without duplicates

    Details:
    A duplicated index means that there are two or more rows for a the
    variant defining that index.  This happens for example with the
    "Overlapped Gene" feature in near_gens.txt annotation file of SNPnexus
    when the variant is in an overlap of multiple genes.  In such cases the
    function collapses the list of gene names into a scalar of "sep" separated
    string of names.

    In general only "object" dtype columns are collapsed into a "sep"
    separated string.  For other dtypes simply the first point of the
    duplicates is used and the rest of the points are discarded.
    '''
    # return annot unchanged unless its index has duplicates
    if not any(annot.index.duplicated()):
        return(annot)
    # get the set of index values (variants)
    A = set(annot.index)
    # do something to the row(s) marked by a variant
    def do_var(var):
        lines = annot.loc[[var]].copy()
        # if it's just a single row return it unchanged
        if len(lines) == 1:
            return(lines)
        # otherwise collapse the multiple rows into a single row
        else:
            line = lines.iloc[[0]].copy()
            for col in annot.columns:
                if lines[col].dtype == 'object':
                    line[col] = sep.join(list(lines[col]))
            return(line)
    l = [do_var(a) for a in A]
    val = pd.concat(l, axis=0)
    return(val)

def get_multi_annotations(sample, annotlist, annotdirpath='/home/attila/projects/bsm/results/2020-09-07-annotations'):
    pass

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
