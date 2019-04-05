#! /usr/bin/env python3

import pandas as pd
import numpy as np
import itertools
import os

def make_genotypes():
    '''
    Make all possible joint genotypes of the 4 CEPH/Utah grandparents

    Returns:
    gt: a list of tuples where each tuple is a joint genotype
    '''
    combinations = list(itertools.combinations_with_replacement((0, 1, 2), 4))
    nested_list = [sorted(list({y: y for y in
        itertools.permutations(c)}.keys())) for c in combinations]
    genotypes = list(itertools.chain(*nested_list))[1:]
    return(genotypes)


def get_mixing_ratios(csv='/home/attila/projects/bsm/tables/ceph-dna-mix.csv'):
    '''
    Read mixing ratios for CEPH/Utah grandparent genomic DNAs from CSV into a
    data frame.
    '''
    mixing_ratios = pd.read_csv(csv, index_col='genome')
    mixing_ratios = mixing_ratios.astype('int64')
    mixing_ratios = mixing_ratios.loc['NA12889':'NA12892'] # remove bottom line total
    return(mixing_ratios)

def get_aaf_of_gt(gt=make_genotypes(), mr=get_mixing_ratios()):
    '''
    Alternative allele frequency (in percent) as a function of
    genotype vector
    Because all values of genotypes are even their halves can be
    represented as int64

    Parameters:
    gt: a list of tuples where each tuple is a joint genotype (see make_genotypes)
    mr: the mixing ratios read by get_mixing_ratios

    Returns:
    the mapping from genotypes to aaf(%)
    '''
    def helper(mix):
        ratios = mr[mix]
        freq = {y: sum(y * ratios) / 2 for y in gt}
        gtstr = [''.join([str(y) for y in x]) for x in gt]
        values = [freq[k] for k in gt]
        df = pd.DataFrame({mix: values}, index=gtstr, dtype='int64')
        return(df)
    aaf_of_gt = pd.concat([helper(m) for m in mr.columns], axis='columns')
    return(aaf_of_gt)

def get_gt_of_aaf(aaf_of_gt=get_aaf_of_gt()):
    '''
    Get mapping from aaf to genotype(s)

    Parameters
    aaf_of_gt: mapping from genotypes to aaf(%)

    Returns:
    mapping from aaf(%) to genotypes
    '''
    def helper(aaf, mix):
        res = list(aaf_of_gt[aaf_of_gt[mix] == aaf].index)
        return(res)
    aaf_values = {m: sorted(list(set(aaf_of_gt[m]))) for m in aaf_of_gt.columns}
    gt_of_aaf = {m: {f: helper(f, m) for f in aaf_values[m]} for m in aaf_of_gt.columns}
    #gt_of_aaf = {m: [{f: helper(f, m)} for f in aaf_values[m]] for m in aaf_of_gt.columns}
    return(gt_of_aaf)

def write_gt_of_aaf(mix='mix1',
        dirpath='/home/attila/projects/bsm/results/2019-03-12-prec-recall-design/',
        gt_of_aaf=get_gt_of_aaf()):
    filepath = dirpath + os.sep + 'gt_of_aaf-' + mix
    '''
    Write aaf to genotype(s) mapping into a file

    Parameters
    mix: 'mix1' or 'mix2' or 'mix3'
    dirpath: the directory into which the file will be created
    gt_of_aaf: genotype to aaf mapping

    Returns:
    the path of the file created
    '''
    with open (filepath, 'w') as f:
        #for entry in gt_of_aaf[mix]:
        for k in gt_of_aaf[mix].keys():
            v = gt_of_aaf[mix][k]
            print(k, end='\t', file=f)
            for gt in v:
                print(gt, end='\t', file=f)
            print('', file=f)
    readmepath = dirpath + os.sep + 'gt_of_aaf.info'
    with open (readmepath, 'w') as f:
        readmestr='Mapping from aaf to genotypes'
        print('Mapping from aaf to genotypes', file=f)
        print('aaf(%): column 1', file=f)
        print('genotypes: column 2,...', file=f)
    return(filepath)
