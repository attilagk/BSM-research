import scipy.stats
import numpy as np
import pandas as pd
import os.path
from bsmcalls import readVCF
from bsmcalls import preprocessing

cmc_clinical_synid = 'syn2279441'
cmc_clinical_path = '/home/attila/projects/bsm/resources/CMC_Human_clinical_metadata.csv'
cmc_ancestry_path='/home/attila/projects/bsm/resources/cmc-ancestry/CMC_MSSM-Penn-Pitt_DNA_GENOTYPE_ANCESTRY_GemTools.tsv'

v1 = ['AF', 'ALT', 'BaseQRankSum', 'DP', 'FILTER/PASS', 'FS', 'GWASpval', 'REF', 'ReadPosRankSum', 'SOR', 'VQSLOD', 'chromatinState_DLPFC', 'culprit', 'evolConstrain', 'szdbCNVcount']
v2 = ['Dx', 'AntipsychAtyp', 'AntipsychTyp', 'Institution', 'EV.3']

def read_clinical(ancestry=True):
    # CMC_Human_clinical_metadata.csv
    if not os.path.exists(cmc_clinical_path):
        import synapseclient
        syn = synapseclient.login()
        wdir = '/home/attila/projects/bsm/resources/'
        clinical_syn = syn.get('syn2279441', downloadLocation=wdir, ifcollision='overwrite.local')
        fpath = clinical_syn.path
    else:
        fpath = cmc_clinical_path
    clinical = pd.read_csv(fpath, index_col='Individual ID')
    if ancestry:
        ancestry = pd.read_csv(cmc_ancestry_path, sep='\t', index_col='Individual_ID')
        ancestry = ancestry.drop(columns=['Genotyping_Sample_ID', 'Cluster'])
        clinical = pd.concat([clinical, ancestry], axis=1)
    return(clinical)

def clin_drop(clin, calls, columns=[]):
    '''
    Drop columns from clinical data frame

    clin: data frame with clinical data
    calls: data frame with variant data
    columns: additional columns to drop
    '''
    bsm_indiv = set(calls.index.get_level_values('Individual ID'))
    cmc_indiv = set(clin.index)
    clin = clin.loc[bsm_indiv, :]
    nobs = clin.count().sort_values()
    min_obs = nobs['EV.1']
    cols2drop = list(nobs[nobs < min_obs].index) + columns
    clin = clin.drop(columns=cols2drop)
    return(clin)

def get_data(merge=False, cleancalls=True, categorize=True, cols2drop=['Sex']):
    '''
    Get all data for BSM project: calls and clinical data

    Arguments
    merge: wether to merge calls and clinical data
    cols2drop: list of columns to drop in clinical data

    Value: calls and clinical data frame separately (in a tuple) or merged
    into a single data frame
    '''
    clin = read_clinical()
    calls = readVCF.readVCFs(clean=cleancalls)
    if categorize:
        calls = preprocessing.convert2categorical(calls)
        clin = preprocessing.convert2categorical(clin)
    clin = clin_drop(clin, calls, columns=cols2drop)
    clin = preprocessing.drop_unused_categories(clin)
    calls = preprocessing.drop_unused_categories(calls)
    data = (calls, clin)
    if merge:
        data = merge_data(calls,clin)
    return(data)

def merge_data(calls, clin):
    calls, clin = calls.align(clin, level='Individual ID', axis=0)
    data = pd.concat([calls, clin], axis=1)
    return(data)

def agg_calls_numeric(calls, funlist=[np.mean, np.std]):
    grouped = calls.select_dtypes(include=['float64', 'int64']).groupby('Individual ID')
    val = grouped.agg(funlist)
    return(val)

def agg_calls_categ(col, calls):
    s = calls[col]
    global_mode = scipy.stats.mode(s).mode[0]
    grouped = s.groupby('Individual ID')
    def Marg_mode(x):
        return(global_mode)
    def Frequency(x):
        val = x.value_counts(normalize=True).loc[global_mode]
        return(val)
    def Entropy(x):
        pk = x.value_counts(normalize=True)
        val = scipy.stats.entropy(pk)
        return(val)
    funlist = [Marg_mode, Frequency, Entropy]
    val = grouped.agg(funlist)
    iterables = [[col], ['Marg_mode', 'Frequency', 'Entropy']]
    val.columns = pd.MultiIndex.from_product(iterables, names=['Variable', 'Transform'])
    return(val)

def agg_calls(calls):
    count = pd.DataFrame({'nCalls': calls.groupby('Individual ID').size()})
    numeric = agg_calls_numeric(calls, funlist=[np.mean, np.std])
    val = pd.concat([count, numeric], axis=1)
    return(val)
