#! /usr/bin/env python3

import pandas as pd
import numpy as np
import itertools

def make_genotypes():
    combinations = list(itertools.combinations_with_replacement((0, 1, 2), 4))
    nested_list = [sorted(list({y: y for y in
        itertools.permutations(c)}.keys())) for c in combinations]
    genotypes = list(itertools.chain(*nested_list))[1:]
    return(genotypes)


def get_mixing_ratios(csv='/home/attila/projects/bsm/tables/ceph-dna-mix.csv'):
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

def get_gt_of_aaf(aaf, mix='mix1', aaf_of_gt=get_aaf_of_gt()):
    res = list(aaf_of_gt[aaf_of_gt[mix] == aaf].index)
    return(res)
