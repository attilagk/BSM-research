import scipy.stats
import numpy as np
import pandas as pd
import os.path
from bsmcalls import readVCF
from bsmcalls import preprocessing

cmc_clinical_synid = 'syn2279441'
cmc_clinical_path = '/home/attila/projects/bsm/resources/CMC_Human_clinical_metadata.csv'
cmc_ancestry_path = '/home/attila/projects/bsm/resources/cmc-ancestry/CMC_MSSM-Penn-Pitt_DNA_GENOTYPE_ANCESTRY_GemTools.tsv'
walsh_gsub_path = '/home/attila/projects/bsm/resources/walsh-manifests/genomics_subject02_template_WalshParkASD-corr.csv'
walsh_vcfs_path = '/home/attila/projects/bsm/results/calls/filtered-vcfs-Walsh.tsv'
chess_vcfs_path = '/home/attila/projects/bsm/results/calls/filtered-vcfs.tsv'

v1 = ['AF', 'ALT', 'BaseQRankSum', 'DP', 'FILTER/PASS', 'FS', 'GWASpval', 'REF', 'ReadPosRankSum', 'SOR', 'VQSLOD', 'chromatinState_DLPFC', 'culprit', 'szdbCNVcount']
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
    #clinical.columns = pd.MultiIndex.from_product([['Clinical'], calls.columns], names=['Source', 'Annotation'])
    return(clinical)

def read_walsh_clinical(clin=read_clinical(), indIDs=np.loadtxt(walsh_vcfs_path, dtype=str)[:, 0]):
    '''
    Read genomic_subject for the Walsh data and format it like CMC clinical
    '''
    walshclin = clin.reindex(indIDs)
    walsh_gsub = pd.read_csv(walsh_gsub_path, index_col='src_subject_id', skiprows=1)
    walsh_gsub = walsh_gsub.loc[indIDs]
    # Institution
    walshclin['Institution'] = walsh_gsub['biorepository']
    # Reported Gender
    g = pd.Categorical(walsh_gsub['gender'], categories=['F', 'M'], ordered=False)
    g = g.rename_categories({'F': 'Female', 'M': 'Male'}, inplace=False)
    walshclin['Reported Gender'] = g
    # Ethnicity
    l = ['Black or African American', 'White', 'American Indian/Alaska Native', 'Asian']
    d = {'Black or African American': 'African-American', 'White': 'Caucasian'}
    e = pd.Categorical(walsh_gsub['race'], categories=l, ordered=False)
    e = e.rename_categories(d, inplace=False)
    walshclin['Ethnicity'] = e
    # ageOfDeath (convert months to years)
    walshclin['ageOfDeath'] = walsh_gsub['interview_age'] / 12
    # Dx
    d = pd.Categorical(walsh_gsub['phenotype'], categories=['Normal', 'Autism'], ordered=False)
    d = d.rename_categories({'Normal': 'Control', 'Autism': 'ASD'}, inplace=False)
    walshclin['Dx'] = d
    return(walshclin)

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

def get_datasets(merge=True, cleancalls=True, categorize=True, cols2drop=['Sex']):
    '''
    Get all datasets for BSM project: calls and clinical data for the Chess and Walsh datasets

    Arguments
    merge: wether to merge calls and clinical data
    cols2drop: list of columns to drop in clinical data

    Value: calls and clinical data frame separately (in a tuple) or merged
    into a single data frame
    '''
    calls_chess = readVCF.readVCFs(vcflistpath=chess_vcfs_path, clean=True)
    calls_walsh = readVCF.readVCFs(vcflistpath=walsh_vcfs_path, clean=True)
    clin_chess = read_clinical()
    clin_chess = clin_drop(clin_chess, calls_chess, columns=cols2drop)
    clin_walsh = read_walsh_clinical(clin=clin_chess)
    clin_chess['Dataset'] = 'Chess'
    clin_walsh['Dataset'] = 'Walsh'
    calls = pd.concat([calls_chess, calls_walsh], axis=0)
    clin = pd.concat([clin_chess, clin_walsh[clin_chess.columns]], axis=0)
    if categorize:
        calls = preprocessing.convert2categorical(calls)
        clin = preprocessing.convert2categorical(clin)
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

def agg_calls_categ(calls):
    def do_col(col, df):
        s = calls[col]
        global_mode = scipy.stats.mode(s).mode[0]
        grouped = s.groupby('Individual ID')
        def marg_mode(x):
            return(global_mode)
        def frequency(x):
            val = x.value_counts(normalize=True).loc[global_mode]
            return(val)
        def entropy(x):
            pk = x.value_counts(normalize=True)
            val = scipy.stats.entropy(pk)
            return(val)
        funlist = [marg_mode, frequency, entropy]
        val = grouped.agg(funlist)
        iterables = [[col], ['marg_mode', 'frequency', 'entropy']]
        val.columns = pd.MultiIndex.from_product(iterables, names=['Variable', 'Transform'])
        return(val)
    df = calls.select_dtypes(include=['category'])
    l = [do_col(c, df) for c in df.columns]
    agg_calls = pd.concat(l, axis=1)
    return(agg_calls)

def agg_calls(calls):
    arrays = [['nCalls'], ['count']]
    tuples = list(zip(*arrays))
    index = pd.MultiIndex.from_tuples(tuples, names=['Variable', 'Transform'])
    count = pd.DataFrame({'nCalls': calls.groupby('Individual ID').size()})
    count.columns = index
    numeric = agg_calls_numeric(calls, funlist=[np.mean, np.std])
    categ = agg_calls_categ(calls)
    val = pd.concat([count, numeric, categ], axis=1)
    return(val)

def get_nsamples(df, margin=False):
    '''
    Get number of samples for each Dx

    Arguments
    df: a calls or annot-like DataFrame
    margin: whether to sum across all Dx to get all individuals

    Value:
    A dict with Dx values as keys and the numbers of individuals as values.
    '''
    gb_list = ['Tissue', 'Dx'] if 'Tissue' in df.index.names else 'Dx'
    d = {name: len(group.index.get_level_values('Individual ID').unique()) for name, group in df.groupby('Dx')}
    if margin:
        d.update({'All': np.sum(list(d.values()))})
    return(d)
