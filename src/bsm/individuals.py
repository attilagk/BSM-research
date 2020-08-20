import pandas as pd
import os.path
from bsm import readVCF

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

def get_data(merge=False, cols2drop=[]):
    '''
    Get all data for BSM project: calls and clinical data

    Arguments
    merge: wether to merge calls and clinical data
    cols2drop: list of columns to drop in clinical data

    Value: calls and clinical data frame separately (in a tuple) or merged
    into a single data frame
    '''
    clin = read_clinical()
    calls = readVCF.readVCFs()
    clin = clin_drop(clin, calls, columns=cols2drop)
    clin = impute(clin)
    data = (calls, clin)
    if merge:
        data = merge_data(calls,clin)
    return(data)

def merge_data(calls,clin):
    calls, clin = calls.align(clin, level='Individual ID', axis=0)
    data = pd.concat([calls, clin], axis=1)
    return(data)
