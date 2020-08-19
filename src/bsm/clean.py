import readVCF
import pandas as pd
import re

v1 = ['AF', 'ALT', 'BaseQRankSum', 'DP', 'FILTER/PASS', 'FS', 'GWASpval', 'REF', 'ReadPosRankSum', 'SOR', 'VQSLOD', 'chromatinState_DLPFC', 'culprit', 'evolConstrain', 'szdbCNVcount']
v2 = ['Dx', 'AntipsychAtyp', 'AntipsychTyp', 'Institution', 'EV.3']

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

