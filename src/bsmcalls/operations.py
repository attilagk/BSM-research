import pandas as pd
import numpy as np
import functools

def series_of_sets_intersect(series_of_sets, query):
    if not isinstance(query, set):
        query = {query}
    s = series_of_sets.dropna()
    s = s.apply(lambda x: len(x.intersection(query)))
    s = s.reindex_like(series_of_sets).fillna(0).astype(np.int64)
    return(s)

def dicts2sets(series_of_dicts, listvalued=True):
    def fun_listvalued(x):
        val = functools.reduce(lambda a, b: a + b, x.values())
        val = set(val)
        return(val)
    def fun_scalarvalued(x):
        val = set(x.values())
        return(val)
    fun = fun_listvalued if listvalued else fun_scalarvalued
    s = series_of_dicts.dropna()
    s = s.apply(fun)
    s = s.reindex_like(series_of_dicts)
    return(s)

features = {'near_gens_Overlapped Gene': lambda x: x,
            'near_gens_Type': lambda x: dicts2sets(x, True),
            'near_gens_Annotation': lambda x: dicts2sets(x, True),
            'tfbs_TFBS Name': lambda x: x,
            'tfbs_TFBS Accession': lambda x: dicts2sets(x, False),
            'tarbase_miRNA': lambda x: x,
            'tarbase_Accession': lambda x: dicts2sets(x, False),
            }

def query(queryitems, feature, data):
    if isinstance(queryitems, list):
        querydict = dict(zip(queryitems, queryitems))
    else:
        querydict = queryitems
    fun = features[feature]
    series_of_sets = fun(data[feature])
    l = [series_of_sets_intersect(series_of_sets, q) for q in querydict.values()]
    df = pd.concat(l, axis=1)
    df.columns = querydict.keys()
    return(df)
