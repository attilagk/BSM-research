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

# TODO: GTEx

annotlist = ['gen_coords', 'ensembl', 'near_gens', # Gene Annotation
                      'sift', 'polyphen', # Protein Effect
                      '1KGen', 'gnomad', # Population Data
                      'tfbs', 'vista', 'cpg', 'targetscan', 'tarbase', 'encode', 'roadmap', 'regbuild', # Regulatory Elements
                      'clinvar', # Phenotype/Disease Association
                      'phast', 'gerp', # Conserved Elements
                      'structvar', # Structural Variations
                      'cadd', 'fitcons', 'eigen', 'fathmm', 'gwava', 'deepsea', 'funseq2', 'remm'] # Non-coding Variation Scoring

annotlists = {'Gene Annotation': ['gen_coords', 'ensembl', 'near_gens'],
              'Protein Effect': ['sift', 'polyphen'],
              'Population Data': ['1KGen', 'gnomad'],
              'Regulatory Elements': ['tfbs', 'vista', 'cpg', 'targetscan', 'tarbase', 'encode', 'roadmap', 'regbuild'],
              'Phenotype/Disease Association': ['clinvar'],
              'Conserved Elements': ['phast', 'gerp'],
              'Structural Variations': ['structvar'],
              'Non-coding Variation Scoring': ['cadd', 'fitcons', 'eigen', 'fathmm', 'gwava', 'deepsea', 'funseq2', 'remm']}

columns2drop = ['ensembl_Variant', 'ensembl_Strand', 'ensembl_CDNA Position', 'ensembl_CDS Position', 'ensembl_Proteins',
                'ensembl_Symbol', 'ensembl_Gene', 'ensembl_Transcript', 'ensembl_Predicted Function', 'ensembl_AA Position',
                'ensembl_AA Change', 'ensembl_Detail', 'ensembl_Splice Distance',
                'sift_dbSNP', 'sift_Variant', 'sift_Transcript',
                'polyphen_dbSNP', 'polyphen_Variant', 'polyphen_Transcript',
                'polyphen_Gene', 'polyphen_AA Position', 'polyphen_Wild AA', 'polyphen_Mutant AA',
                '1KGen_dbSNP', '1KGen_REF Allele', '1KGen_ALT Allele', '1KGen_Minor Allele',
                'tarbase_Gene', 
                'encode_Region Start', 'encode_Region End', 'encode_Feature Type Class', 'encode_Feature Type', 'encode_Epigenome',
                'roadmap_Region Start', 'roadmap_Region End', 'roadmap_Feature Type Class', 'roadmap_Feature Type', 'roadmap_Epigenome',
                'regbuild_Region Start', 'regbuild_Region End', 'regbuild_Feature Type', 'regbuild_Epigenome', 'regbuild_Activity',
                'clinvar_Variation', 'clinvar_Type',
                'phast_Region Start', 'phast_Region End', 'phast_Id',
                'gerp_Region Start', 'gerp_Region End',
                'structvar_Reference', 'structvar_Chrom Start', 'structvar_Chrom End', 'structvar_Type', 'structvar_Method', 'structvar_Sample', 'structvar_Gain',
                'clinvar_Type',
                'cadd_Variant', 'fitcons_Region Start', 'fitcons_Region End',
                'eigen_Variant', 'fathmm_Variant', 'gwava_Known SNP','gwava_Region Score',  'gwava_TSS Score', 'gwava_Unmatched Score',
                'deepsea_Variant', 'funseq2_Variant'
                ]

columns2float = ['gen_coords_Minor Allele Global Frequency',
                 'sift_Score', 'polyphen_Score',
                 '1KGen_AFR Frequency', '1KGen_AMR Frequency', '1KGen_EAS Frequency', '1KGen_EUR Frequency', '1KGen_SAS Frequency',
                 'gnomad_AFR Frequency', 'gnomad_AMR Frequency', 'gnomad_ASJ Frequency',
                 'gnomad_EAS Frequency', 'gnomad_FIN Frequency', 'gnomad_NFE Frequency',
                 'gnomad_OTH Frequency', 'gnomad_SAS Frequency',
                 'cpg_CpG %', 'cpg_C/G %', 'cpg_Ratio',
                 'gerp_Element RS Score', 'gerp_Base RS Score',
                 'cadd_Raw Score', 'cadd_PHRED', 'fitcons_Fitness Score',
                 'eigen_Score', 'eigen_PC Score',
                 'fathmm_Non-coding Score', 'fathmm_Coding Score',
                 'deepsea_Functional Significance Score', 'deepsea_eQTL Probability',
                 'deepsea_GWAS Probability', 'deepsea_HGMD Probability',
                 'funseq2_Non-coding Score', 'remm_ReMM Score'
                   ]

columns2integer = ['near_gens_Distance to Nearest Upstream Gene', 'near_gens_Distance to Nearest Downstream Gene',
                   'gen_coords_Contig Position',
                   'sift_AA Position',
                   'tfbs_Region Start', 'tfbs_Region End',
                   'vista_Region Start', 'vista_Region End', 'vista_Score',
                   'cpg_Region Start', 'cpg_Region End', 'cpg_Length',
                   'targetscan_Region Start', 'targetscan_Region End', 'targetscan_Score', 'targetscan_Strand',
                   'tarbase_Region Start', 'tarbase_Region End', 'tarbase_Strand',
                   'phast_Score', 'structvar_PubMed', 'structvar_Loss',
                   ]

columns2split_keep1st = ['sift_Gene', 'sift_Wild AA', 'sift_Mutant AA', 'sift_Prediction',
                         'polyphen_Prediction',
                         ]


aa_alphabet = list('ACDEFGHIKLMNPRSTVWY')

columns2categorize = {'gen_coords_dbSNP': None,
                      'gen_coords_REF Allele': list('ACGT'),
                      'gen_coords_ALT Allele (IUPAC)': list('ACGT'),
                      'gen_coords_Minor Allele': list('ACGT'),
                      'gen_coords_Contig': None,
                      'gen_coords_Band': None,
                      'near_gens_Nearest Upstream Gene': None,
                      'near_gens_Type of Nearest Upstream Gene': None,
                      'near_gens_Nearest Downstream Gene': None,
                      'near_gens_Type of Nearest Downstream Gene': None,
                      'sift_Gene': None,
                      'sift_Wild AA': aa_alphabet,
                      'sift_Mutant AA': aa_alphabet,
                      'sift_Prediction': None,
                      'polyphen_Prediction': None,
                      'gnomad_dbSNP': None,
                      'gnomad_REF Allele': list('ACGT'),
                      'gnomad_ALT Allele': list('ACGT'),
                      'gnomad_Minor Allele': list('ACGT'),
                      'clinvar_Clinical Significance': None,
                      'clinvar_Phenotypes': None,
                      'vista_Vista Item': None,
                      'fathmm_Non-coding Groups': ['AB', 'A', 'AC', 'ALL', 'ADB', 'AD', 'ABC', 'ADC', 'D', 'B'],
                      'fathmm_Coding Groups': None
                      }

na_values = dict()

none_colon_none='^None(:None)*$'

NA2remove = {'near_gens_Overlapped Gene': ['None'],
             'near_gens_Type': ['None'],
             'near_gens_Annotation': ['None'],
             'near_gens_Nearest Upstream Gene': none_colon_none,
             'near_gens_Type of Nearest Upstream Gene': none_colon_none,
             'near_gens_Nearest Downstream Gene': none_colon_none,
             'near_gens_Type of Nearest Downstream Gene': none_colon_none,
             }


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
                          vcflistpath='/big/results/bsm/calls/filtered-vcfs-Chess-Walsh.tsv',
                          annotdirpath='/home/attila/projects/bsm/results/2020-09-07-annotations',
                          na_values={}, simplecolumns=True):
    vcflist = pd.read_csv(vcflistpath, sep='\t', names=['sample', 'file'], index_col='sample')
    samplestr = '((MSSM|PITT)_[0-9]+)_(NeuN_pl|NeuN_mn|muscle)'
    def sample2indivID(sample):
        return(re.sub(samplestr, 'CMC_\\1', sample))
    def sample2tissue(sample):
        #if re.match('.*Walsh.*', vcflistpath):
        if not re.match(samplestr, sample):
            return('frontal cortex') # Walsh data
        return(re.sub(samplestr, '\\3', sample)) # Chess data
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
        try:
            annot = pd.concat([get_annot(s, annotyp) for s in vcflist.index], axis=0)
        except ValueError:
            annot = None
        return(annot)
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


def create_colsdict():
    '''
    Create input dictionary for regularize_categ_cols
    '''
    colsdict = {}
    # order reflecting severity of effect
    l = ['Deleterious', 'Deleterious - Low Confidence', 'Tolerated', 'Tolerated - Low Confidence']
    colsdict.update({'sift_Prediction': l})
    # order reflecting increasing frequency of categories in the data set
    l = ['Polymerase', 'Open Chromatin', 'Transcription Factor', 'Histone']
    colsdict.update({'encode_Feature Type Class': l})
    l = ['intronic (splice_site)', 'coding', 'intronic', '5utr', '3utr', '5upstream', '3downstream', 'non-coding intronic', 'non-coding']
    colsdict.update({'ensembl_Predicted Function': l})
    def read_categories(fpath):
        with open(fpath) as f:
            val = f.readlines()
            val = [x.strip() for x in val]
            return(val)
    regbuild_epigenomes = read_categories('/big/results/bsm/2020-09-07-annotations/regbuild-epigenomes')
    colsdict.update({'regbuild_Epigenome': regbuild_epigenomes})
    colsdict.update({'structvar_Type': ['complex', 'loss', 'gain']})
    return(colsdict)


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
    '''
    Obsoleteed main function superseded by merge_snpnexus_with_other_annotations
    '''
    annot = get_multi_annotations(annotlist)
    numeric_cols = annot.select_dtypes(np.number).columns
    cols2binarize = [c for c in numeric_cols if c in cols2process]
    annot = binarize_cols(cols2binarize, annot, calls, suffix='_bin')
    return(annot)


def load_data(picklepath='/home/attila/projects/bsm/results/2020-09-07-annotations/annotated-calls.p'):
    '''
    Load annotated calls from pickle file
    '''
    with open(picklepath, 'rb') as f:
        data = pickle.load(f)
    return(data)


def filter_fulldata(fulldata):
    '''
    Filter all VCF records: keep only NeuN_pl from Chess data and HC from Walsh data
    '''
    chess = fulldata.loc[fulldata.Dataset == 'Chess'].xs('NeuN_pl', level='Tissue')
    HC_list = ['HC/PASS', 'HC;PASS/PASS']
    walsh = fulldata.loc[fulldata.Dataset == 'Walsh'].xs('frontal cortex', level='Tissue')
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

def postprocess_annot(annot, cols2drop=columns2drop, cols2float=columns2float,
                      columns2split_keep1st=columns2split_keep1st, 
                      cols2integer=columns2integer, NA2remove=NA2remove,
                      cols2categorize=columns2categorize):
    # drop unnecessary columns
    try:
        annot = annot.drop(labels=cols2drop, axis=1)
    except KeyError:
        annot = annot.copy()
    for c in cols2float + cols2integer:
        annot[c] = pd.to_numeric(annot[c], errors='coerce')
    for c in cols2integer:
        annot[c] = pd.Series(annot[c], dtype=pd.Int64Dtype())
    for k, v in NA2remove.items():
        if isinstance(v, list):
            b = annot[k].isin(v)
        else:
            b = annot[k].str.match(v)
        annot.loc[b, k] = np.nan
    for c in columns2split_keep1st:
        annot[c] = annot[c].str.extract('^([^:]+):.*', expand=False)
    for k, v in columns2categorize.items():
        if v is None:
            s = pd.Categorical(annot[k])
        else:
            s = pd.Categorical(annot[k], categories=v)
        annot[k] = s
    # odds and ends
    annot['fitcons_p-Value'] = annot['fitcons_p-Value'].str.strip('<').astype('float64')
    s = annot['cpg_CpG Island'].str.strip('CpG: ')
    s = pd.to_numeric(s, errors='coerce')
    annot['cpg_CpG Island'] = pd.Series(s, dtype=pd.Int64Dtype())
    return(annot)

def merge_snpnexus_with_other_annotations(annotlist=annotlist,
                                          na_values=na_values,
                                          colsdict=create_colsdict(),
                                          fpath='/home/attila/projects/bsm/results/2020-09-07-annotations/annot.p',
                                          calls=individuals.get_datasets()):
    '''
    Main function: read SNPnexus annotations for the full Chess and Walsh datasets and merge them with other annotations
    '''
    if os.path.exists(fpath):
        print('loading annot DataFrame from', fpath)
        with open(fpath, 'rb') as f:
            annot = pickle.load(f)
    else:
        vcflistpath = '/home/attila/projects/bsm/results/calls/filtered-vcfs-Chess-Walsh.tsv'
        annotdirpath = '/home/attila/projects/bsm/results/2020-09-07-annotations'
        annot = get_multi_annotations(annotlist, vcflistpath, annotdirpath, na_values)
        pickle.dump(annot, open(fpath, 'wb'))
    return(annot)
    
    annot = regularize_categ_cols(colsdict, annot, calls, nafillval='other')
    #s = annot['regbuild_Epigenome']
    #annot['regbuild_Epigenome_nervoussys_bin'] = np.int8(s.isin(regbuild_epigenomes[:7]))
    data = pd.concat([calls, annot], axis=1)
    return(data)
