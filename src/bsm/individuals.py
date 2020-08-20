import pandas as pd
import re
import os.path

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

def impute(clin):
    nobs = clin.count()
    max_obs = len(clin)
    cols2impute = set(nobs[nobs < max_obs].index)
    # drop all non-numeric columns with missing data
    numeric = set(clin.select_dtypes(include=['float64', 'int64']).columns)
    cols2drop = cols2impute - cols2impute.intersection(numeric)
    clin = clin.drop(columns=cols2drop)
    # function to calculate impute value
    def helper(col):
        return(clin[col].mean())
    # do the imputation for all columns
    for col in cols2impute:
        clin[col].fillna(value=helper(col), inplace=True)
    return(clin)

def preselect(calls, vnames=v1 + v2):
    calls = calls.loc[:, vnames]
    return(calls)

def prettify_colnames(calls, repl='_', pattern='[ ./\:	]+'):
    calls = calls.rename(lambda y: re.sub(pattern, repl, y), axis='columns')
    return(calls)

def dummify_var(calls, vname='Dx'):
    s = calls[vname].astype('category')
    if len(s.cat.categories) != 2:
        raise TypeError(vname + ' is not a binary variable')
    calls[vname] = s.cat.rename_categories([0, 1]).astype('int32')
    return(calls)

def impute_vars(calls, vnames=['ReadPosRankSum', 'EV_3'], v1=v1, v2=v2):
    repl = '_'
    pattern = '[ ./\:	]+'
    v1 = [re.sub(pattern, repl, y) for y in v1]
    v2 = [re.sub(pattern, repl, y) for y in v2]
    def helper(vname):
        if vname in set(v1):
            val = calls.mean()
        elif vname in set(v2):
            val = calls.groupby('Individual ID')[vname].first().mean()
        else:
            raise ValueError("Couldn't find", vname)
        return(val)
    for v in vnames:
        print(v)
        calls[v].fillna(value=helper(v), inplace=True)
    return(calls)


def standardize(calls):
    stdcalls = calls.apply(lambda y: (y - y.mean()) / y.std() if (y.dtype == 'float64' or y.dtype == 'int64') else y, axis=0)
    return(stdcalls)

