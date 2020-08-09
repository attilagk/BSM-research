import readVCF
import pandas as pd
import statsmodels.api as sm

selvars = ['AntipsychAtyp', 'Alcohol', 'AntipsychTyp', 'Ethnicity']
cmc_clinical_path = '/big/results/bsm/2020-08-05-cmc-clinical/CMC_Human_clinical_metadata.csv'
ancestry_path='/home/attila/projects/bsm/resources/cmc-ancestry/CMC_MSSM-Penn-Pitt_DNA_GENOTYPE_ANCESTRY_GemTools.tsv'

def clean_clinical(calls, remove_ancestry=False):
    clin = pd.read_csv(cmc_clinical_path)
    todrop = set(calls.columns).intersection(set(clin.columns)) - set(selvars)
    if remove_ancestry:
        ancestry = pd.read_csv(ancestry_path, sep='\t', index_col='Individual_ID')
        todrop = todrop.union(set(ancestry.columns))
    calls = calls.drop(columns=todrop)
    return(calls)
