import statsmodels as sm
import re


def prettify_colnames(calls, repl='_', pattern='[ ./\:	]+'):
    calls = calls.rename(lambda y: re.sub(pattern, repl, y), axis='columns')
    return(calls)

def dummify_var(calls, vname='Dx'):
    s = calls[vname].astype('category')
    if len(s.cat.categories) != 2:
        raise TypeError(vname + ' is not a binary variable')
    calls[vname] = s.cat.rename_categories([0, 1]).astype('int32')
    return(calls)

def standardize(calls):
    stdcalls = calls.apply(lambda y: (y - y.mean()) / y.std() if (y.dtype == 'float64' or y.dtype == 'int64') else y, axis=0)
    return(stdcalls)

