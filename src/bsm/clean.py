import readVCF
import pandas as pd
import re

cmc_clinical_synid = 'syn2279441'
cmc_clinical_path = '/home/attila/projects/bsm/resources/CMC_Human_clinical_metadata.csv'
cmc_ancestry_path='/home/attila/projects/bsm/resources/cmc-ancestry/CMC_MSSM-Penn-Pitt_DNA_GENOTYPE_ANCESTRY_GemTools.tsv'

v1 = ['AF', 'ALT', 'BaseQRankSum', 'DP', 'FILTER/PASS', 'FS', 'GWASpval', 'REF', 'ReadPosRankSum', 'SOR', 'VQSLOD', 'chromatinState_DLPFC', 'culprit', 'evolConstrain', 'szdbCNVcount']
v2 = ['Dx', 'AntipsychAtyp', 'AntipsychTyp', 'Institution', 'EV.3']

def read_clinical():
    # CMC_Human_clinical_metadata.csv
    syn = synapseclient.login()
    wdir = '/home/attila/projects/bsm/resources/'
    clinical_syn = syn.get('syn2279441', downloadLocation=wdir, ifcollision='overwrite.local')
    clinical = pd.read_csv(clinical_syn.path, index_col='Individual ID')

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

def autofilter(calls):
    ourindivs = set(calls.index.get_level_values(0))
    clinical = pd.read_csv(cmc_clinical_path, index_col='Individual ID')
    ancestry = pd.read_csv(cmc_ancestry_path, sep='\t', index_col='Individual_ID')
    set(clinical.index)

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

