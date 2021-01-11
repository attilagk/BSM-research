import pandas as pd
import numpy as np
import re
import os.path
import functools
import operator
import copy
import itertools
import ensembl_rest
import pickle
from bsmcalls import individuals


def tsvpath2annotname(tsvpath):
    val = re.sub('.txt', '', os.path.basename(tsvpath))
    return(val)

def read_TXT_per_annotation(tsvpath, indivID, tissue='NeuN_pl',
                            annotname=None, na_values=[], simplecolumns=True):
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
    annot = pd.read_csv(tsvpath, sep='\t', na_values=na_values)
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
    annot = annot.loc[:, ~annot.columns.isin(['Variation ID', 'Chromosome', 'Position'])]
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
                          na_values={}, simplecolumns=True):
    vcflist = pd.read_csv(vcflistpath, sep='\t', names=['sample', 'file'], index_col='sample')
    samplestr = '((MSSM|PITT)_[0-9]+)_(NeuN_pl|NeuN_mn|muscle)'
    def sample2indivID(sample):
        return(re.sub(samplestr, 'CMC_\\1', sample))
    def sample2tissue(sample):
        if re.match('.*Walsh.*', vcflistpath):
            return('frontal cortex')
        return(re.sub(samplestr, '\\3', sample))
    def get_annot(sample, annotyp):
        sampledir = annotdirpath + os.path.sep + sample
        tsvpath = sampledir + os.path.sep + annotyp + '.txt'
        indivID = sample2indivID(sample)
        tissue = sample2tissue(sample)
        na_val = na_values[annotyp] if annotyp in na_values.keys() else []
        try:
            annot = read_TXT_per_annotation(tsvpath, indivID, tissue,
                                            simplecolumns=simplecolumns, na_values=na_val)
            annot = annotation_duplicates(annot, sep=':')
        except ValueError:
            annot = None
        return(annot)
    def do_annotyp(annotyp):
        a = pd.concat([get_annot(s, annotyp) for s in vcflist.index], axis=0)
        return(a)
    annot = pd.concat([do_annotyp(a) for a in annotlist], axis=1)
    return(annot)

def extended_columns(columns, cols2insert, suffix='_bin'):
    def helper(c):
        val = [c, c + suffix] if c in cols2insert else [c]
        return(val)
    l = [helper(c) for c in columns]
    extcolumns = functools.reduce(operator.concat, l)
    return(extcolumns)

def binarize_cols(cols2binarize, annot, calls, suffix='_bin', do_categ=False):
    '''
    Binarize the selected columns of annot and reindex it to match calls

    Arguments
    cols2binarize: list of column names to binarize
    annot: the pandas DataFrame containing cols2binarize
    calls: annot will be reindexed according to this DataFrame
    suffix: of the names of binarized columns

    Value: a copy of annot extended with the binarized columns
    '''
    val = annot.copy()
    columns = extended_columns(columns=annot.columns, cols2insert=cols2binarize, suffix=suffix)
    val = val.reindex(columns=columns)
    val = val.reindex(index=calls.index)
    def do_col(col):
        s = np.int8([not y for y in val[col].isna()])
        val[col + suffix] = pd.Categorical(s) if do_categ else s
    for col in cols2binarize:
        do_col(col)
    return(val)

def regularize_categ_cols(colsdict, annot, calls, nafillval='other'):
    '''
    Regularize categorical columns in annot: map vectors to scalars and fill NAs

    Arguments
    colsdict: a dictonary of column names (keys) and the list of their ordered categories
    annot: the DataFrame to be modified (copy)
    calls: the DataFrame based on our VCF
    nafillval: the value to replace missing values with

    Value: the modified copy of annot

    Details:
    When there are multiple values for a row/variant in a given column
    represented in colsdict then the corresponding order of categories
    determines which value is kept and which are removed.
    '''
    val = annot.copy()
    # deep copy is necessary to prevent mutation of colsdict
    d = copy.deepcopy(colsdict)
    for col, categories in d.items():
        categories += ['other']
        s = val[col]
        def helper(old):
            if old is np.nan:
                return(old)
            if not isinstance(old, str):
                raise ValueError('expected either str or np.nan')
            l = old.split(':')
            for cat in categories:
                if cat in l: return(cat)
            return(old)
        s = [helper(x) for x in s]
        s = pd.Categorical(s, categories=categories, ordered=True)
        s = s.fillna(nafillval)
        val[col] = s
    return(val)

def do_annot(annotlist, calls, cols2process=None):
    annot = get_multi_annotations(annotlist)
    numeric_cols = annot.select_dtypes(np.number).columns
    cols2binarize = [c for c in numeric_cols if c in cols2process]
    annot = binarize_cols(cols2binarize, annot, calls, suffix='_bin')
    return(annot)

def load_data(picklepath='/home/attila/projects/bsm/results/2020-09-07-annotations/annotated-calls.p'):
    with open(picklepath, 'rb') as f:
        data = pickle.load(f)
    return(data)

def filter_fulldata(fulldata):
    chess = fulldata.loc[fulldata.Dataset == 'Chess'].xs('NeuN_pl', level='Tissue')
    HC_list = ['HC/PASS', 'HC;PASS/PASS']
    walsh = fulldata.loc[fulldata.Dataset == 'Walsh']
    walsh = walsh.loc[walsh['FILTER/PASS'].isin(HC_list)]
    data = pd.concat([chess, walsh], axis=0)
    return(data)

def str2num(annot, colname):
    s = annot[colname].copy()
    s = pd.to_numeric(s, errors='coerce')
    return(s)

def multi_str2num(annot, colnames):
    '''
    The multi feature version of str2list
    '''
    df = annot.copy()
    l = [str2num(df, c) for c in colnames]
    d = dict(zip(colnames, l))
    df = df.assign(**d)
    return(df)

def str2list(annot, colname, nonestr='None', sepstr=':'):
    '''
    Convert a string to a list removing nonestr and splitting on sepstr

    Arguments
    annot: pandas DataFrame returned by merge_snpnexus_with_other_annotations, do_annot or read by load_data
    colname: the column name of the feature
    nonestr: this means NA or similar
    sepstr: splitting should be done on this string

    Value: the engineered copy of annot[colname]; a pd.Series object

    Details:
    Typically used with near_gens_Overlapped Gene, near_gens_Type,
    near_gens_Annotation as well as the upstream and downstream version of
    these features.

    '''
    s = annot[colname].copy()
    s = s.fillna('')
    s = s.str.replace(',', sepstr)
    pattern = '^None(' + sepstr + 'None)*$'
    s = s.apply(lambda x: re.sub(pattern, '', x))
    s = s.str.split(sepstr)
    return(s)


def multi_str2list(annot, colnames, nonestr='None', sepstr=':'):
    '''
    The multi feature version of str2list

    Arguments are similar to those of str2list except...
    colnames: the list of column name of the feature s

    Value: the engineered copy of annot; a pd.DataFrame object
    '''
    df = annot.copy()
    l = [str2list(df, c, nonestr=nonestr, sepstr=sepstr) for c in colnames]
    d = dict(zip(colnames, l))
    df = df.assign(**d)
    return(df)

def ensembl_description(annot, colname='near_gens_Overlapped Gene set'):
    def helper(symbolset=annot[colname][1]):
        if len(symbolset) == 0:
            return(dict())
        d = dict()
        for sym in symbolset:
            try:
                description = ensembl_rest.symbol_lookup('homo sapiens', sym, expand=1)['description']
            except Exception:
                description = 'No valid lookup found for symbol ' + sym
            d[sym] = description
        return(d)
    val = [helper(y) for y in annot[colname]]
    return(val)

def insert_col(s, df, oldname, newname, inplace=False):
    '''
    Insert vector s into df as column newname after column oldname
    '''
    if inplace:
        D = df
    else:
        D = df.copy()
    if newname in D.columns:
        return(D)
    ix = list(D.columns).index(oldname)
    D.insert(ix + 1, column=newname, value=s)
    return(D)

def merge_snpnexus_with_other_annotations(calls=individuals.get_datasets(), testmode=False):
    '''
    Main function: read SNPnexus annotations for the full Chess and Walsh datasets and merge them with other annotations
    '''
    # SNPnexus annotations
    annotlist = ['1KGen', 'cpg', 'deepsea', 'encode', 'ensembl', 'gerp', 'near_gens', 'phast', 'regbuild', 'sift',  'structvar', 'tfbs']
    # custom values
    na_values = {}
    na_values.update({'1KGen': {'AFR Frequency': 'None', 'AMR Frequency': 'None', 'EAS Frequency': 'None', 'EUR Frequency': 'None', 'SAS Frequency': 'None'}})
    # to binarize columns
    cols2binarize = []
    cols2binarize += ['1KGen_AFR Frequency', '1KGen_AMR Frequency', '1KGen_EAS Frequency', '1KGen_EUR Frequency', '1KGen_SAS Frequency']
    cols2binarize += ['cpg_CpG Island']
    cols2binarize += ['gerp_Element RS Score']
    cols2binarize += ['phast_Score']
    cols2binarize += ['tfbs_TFBS Name']
    cols2binarize += ['structvar_Type']
    # to regularize columns
    colsdict = {}
    # order reflecting severity of effect
    l = ['Deleterious', 'Deleterious - Low Confidence', 'Tolerated', 'Tolerated - Low Confidence']
    colsdict.update({'sift_Prediction': l})
    # order reflecting increasing frequency of categories in the data set
    l = ['Polymerase', 'Open Chromatin', 'Transcription Factor', 'Histone']
    colsdict.update({'encode_Feature Type Class': l})
    l = ['coding', 'intronic', 'intronic (splice_site)', '5utr', '3utr', '5upstream', '3downstream', 'non-coding intronic', 'non-coding']
    colsdict.update({'ensembl_Predicted Function': l})
    def read_categories(fpath):
        with open(fpath) as f:
            val = f.readlines()
            val = [x.strip() for x in val]
            return(val)
    regbuild_epigenomes = read_categories('/big/results/bsm/2020-09-07-annotations/regbuild-epigenomes')
    colsdict.update({'regbuild_Epigenome': regbuild_epigenomes})
    colsdict.update({'structvar_Type': ['complex', 'loss', 'gain']})
    if testmode:
        annotlist = annotlist[:2] + ['sift' + 'regbuild']
        cols2binarize = cols2binarize[:6]
        features = ['sift_Prediction', 'regbuild_Epigenome']
        colsdict = {k: colsdict[k] for k in features}
    vcflistpaths = ['/home/attila/projects/bsm/results/calls/filtered-vcfs' + s + '.tsv' for s in ['', '-Walsh']]
    annotdirpath = '/home/attila/projects/bsm/results/2020-09-07-annotations'
    lannot = [get_multi_annotations(annotlist, p, annotdirpath, na_values) for p in vcflistpaths]
    annot = pd.concat(lannot, axis=0)
    annot = binarize_cols(cols2binarize, annot, calls, suffix='_bin')
    # read epigenome names for Ensemble Regulatory Build
    # https://useast.ensembl.org/info/genome/funcgen/regulatory_build.html
    annot = regularize_categ_cols(colsdict, annot, calls, nafillval='other')
    s = annot['regbuild_Epigenome']
    annot['regbuild_Epigenome_nervoussys_bin'] = np.int8(s.isin(regbuild_epigenomes[:7]))
    data = pd.concat([calls, annot], axis=1)
    return(data)
