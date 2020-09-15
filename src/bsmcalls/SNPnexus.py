import pandas as pd
import numpy as np
import re
import os.path
import functools
import operator


def tsvpath2annotname(tsvpath):
    val = re.sub('.txt', '', os.path.basename(tsvpath))
    return(val)

def read_TXT_per_annotation(tsvpath, indivID, tissue='NeuN_pl',
                            annotname=None, simplecolumns=True):
    '''
    Reads a TXT_per_annotation file of SNPnexus into an indexed DataFrame

    Arguments
    tsvpath: path to the file
    indivID: Individual ID without the CMC_ prefix
    tissue: NeuN_pl|NeuN_mn|muscle
    annotname: the name of annotation; if None it is the basename of tsvpath without the .txt extention
    simplecolumns: True if we want to avoid multilevel type column names (MultiIndex)

    Value: an annot DataFrame
    '''
    if annotname is None:
        annotname = tsvpath2annotname(tsvpath)
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
    simpleindex = [annotname + '_' + a for a in annot.columns]
    multiindex = pd.MultiIndex.from_product([[annotname], annot.columns], names=['Source', 'Annotation'])
    annot.columns = simpleindex if simplecolumns else multiindex
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

def get_multi_annotations(annotlist,
                          vcflistpath='/big/results/bsm/calls/filtered-vcfs.tsv',
                          annotdirpath='/home/attila/projects/bsm/results/2020-09-07-annotations',
                          simplecolumns=True):
    vcflist = pd.read_csv(vcflistpath, sep='\t', names=['sample', 'file'], index_col='sample')
    samplestr = '((MSSM|PITT)_[0-9]+)_(NeuN_pl|NeuN_mn|muscle)'
    def sample2indivID(sample):
        return(re.sub(samplestr, 'CMC_\\1', sample))
    def sample2tissue(sample):
        return(re.sub(samplestr, '\\3', sample))
    def get_annot(sample, annotyp):
        sampledir = annotdirpath + os.path.sep + sample
        tsvpath = sampledir + os.path.sep + annotyp + '.txt'
        indivID = sample2indivID(sample)
        tissue = sample2tissue(sample)
        try:
            annot = read_TXT_per_annotation(tsvpath, indivID, tissue, simplecolumns=simplecolumns)
            annot = annotation_duplicates(annot, sep=':')
        except ValueError:
            annot = None
        return(annot)
    def do_annotyp(annotyp):
        a = pd.concat([get_annot(s, annotyp) for s in vcflist.index], axis=0)
        return(a)
    annot = pd.concat([do_annotyp(a) for a in annotlist], axis=1)
    return(annot)

def binarize_cols(cols, annot, calls, suffix='_bin'):
    '''
    Binarize the selected columns of annot and reindex it to match calls

    Arguments
    cols: list of column names to binarize
    annot: the pandas DataFrame containing cols
    calls: annot will be reindexed according to this DataFrame
    suffix: of the names of binarized columns

    Value: a copy of annot extended with the binarized columns
    '''
    val = annot.copy()
    def helper(c):
        val = [c, c + suffix] if c in cols else [c]
        return(val)
    l = [helper(c) for c in annot.columns]
    columns = functools.reduce(operator.concat, l)
    val = val.reindex(columns=columns)
    val = val.reindex(index=calls.index)
    def do_col(col):
        s = np.int8(np.isnan(val[col]))
        val[col + suffix] = pd.Categorical(s)
    for col in cols:
        do_col(col)
    return(val)

def do_annot(annotlist, calls, cols2process=None):
    annot = get_multi_annotations(annotlist)
    numeric_cols = annot.select_dtypes(np.number).columns
    cols2binarize = [c for c in numeric_cols if c in cols2process]
    annot = binarize_cols(cols2binarize, annot, calls, suffix='_bin')
    return(annot)

'''
========================================================================
Specific Annotations
========================================================================
'''
