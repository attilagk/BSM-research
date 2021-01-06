import pandas as pd
import functools
from bsmcalls import SNPnexus

functional_genic = ['coding nonsyn', 'stop-gain', 'intronic (splice_site)']

def read_expand_calls():
    calls = SNPnexus.read_annot()
    colnamel = ['near_gens_Type', 'near_gens_Annotation', 'ensembl_Predicted Function', 'sift_Prediction']
    nonestrl = ['None', 'None', 'other', 'other']
    calls = SNPnexus.expand_multiple_setvalued(calls, colnamel=colnamel, nonestrl=nonestrl)
    return(calls)

def filter_calls(calls, filtl=functional_genic):
    l = [calls[y] for y in filtl]
    boolix = functools.reduce(lambda x, y: x | y, l)
    calls = calls.loc[boolix]
    return(calls)

def get_coding(calls, filtl=None):
    if filtl is None:
        calls = filter_calls(calls, filtl=functional_genic)
    calls = calls.loc[~ (calls['near_gens_Overlapped Gene'] == 'None')]
    return(calls)

def genic_recurrence_df(calls, genesl=None):
    if genesl is None:
        genesl = list(set(calls['near_gens_Overlapped Gene']))
        genesl.sort()
    s = calls.groupby('Individual ID')['near_gens_Overlapped Gene'].agg(set)
    d = {gene: [ind for ind in s.index if gene in s[ind]] for gene in genesl}
    recurrence_ser = pd.Series(d)
    def helper(l):
        val = [sum([calls.xs(key=ind, axis=0, level='Individual ID')['Dx'].unique()[0] == dx for ind in l]) for dx in ['Control', 'SCZ']]
        return(val)
    ss = recurrence_ser.apply(helper)
    recurrence = ss.apply(pd.Series).rename(columns={0: 'n Control indiv', 1: 'n SCZ indiv'})
    recurrence['Individual IDs'] = recurrence_ser
    recurrence = recurrence.sort_values('n SCZ indiv', ascending=False)
    return(recurrence)

def genic_recurrence(calls, filtl=functional_genic, genesl=None):
    coding = get_coding(calls, filtl=filtl)
    recurrence = genic_recurrence_df(coding, genesl=genesl)
    return(recurrence)

def contingency_tab(recurrence):
    l = ['n Control indiv', 'n SCZ indiv']
    contingency = recurrence.value_counts(l).unstack().fillna(0).astype('int64').rename_axis(index='m Control indiv')
    return(contingency)

