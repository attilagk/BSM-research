import pandas as pd
import numpy as np
import functools

_level_names = ['Feature', 'Query']

def series_of_sets_intersect(series_of_sets, query):
    if not isinstance(query, set):
        query = {query}
    s = series_of_sets.dropna()
    s = s.apply(lambda x: len(x.intersection(query)))
    s = s.reindex_like(series_of_sets).fillna(0).astype(np.int8)
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

def scalar2sets(series_of_scalars):
    s = series_of_scalars.dropna().astype(pd.StringDtype())
    s = s.apply(lambda x: set([x]))
    s = s.reindex_like(series_of_scalars)
    return(s)

features = {'near_gens_Overlapped Gene': lambda x: x,
            'near_gens_Type': lambda x: dicts2sets(x, True),
            'near_gens_Annotation': lambda x: dicts2sets(x, True),
            'tfbs_TFBS Name': lambda x: x,
            'tfbs_TFBS Accession': lambda x: dicts2sets(x, False),
            'tarbase_miRNA': lambda x: x,
            'tarbase_Accession': lambda x: dicts2sets(x, False),
            }

def anyquery(feature, data):
    s = data[feature].fillna(0).astype('bool')
    ix = pd.MultiIndex.from_tuples([(feature, 'any')], names=_level_names)
    df = pd.DataFrame(s.to_numpy(), index=s.index, columns=ix)
    return(df)

def query(queryitems, feature, data):
    names = _level_names
    if queryitems is None:
        return(anyquery(feature, data))
    if isinstance(queryitems, list):
        querydict = dict(zip(queryitems, queryitems))
    else:
        querydict = queryitems
    if feature in features.keys():
        fun = features[feature]
    else:
        fun = scalar2sets
    series_of_sets = fun(data[feature])
    l = [series_of_sets_intersect(series_of_sets, q) for q in querydict.values()]
    df = pd.concat(l, axis=1, keys=querydict.keys())
    iterables = [[feature], queryitems]
    ix = pd.MultiIndex.from_product(iterables, names=names)
    df = pd.DataFrame(df.to_numpy(), index=df.index, columns=ix)
    return(df)

def multiquery(querydict, data, do_sum=False, do_sort=False, margin=True):
    l = [query(v, k, data) for k, v in querydict.items()]
    df = pd.concat(l, axis=1)
    if do_sort:
        df = df.sort_index(axis=1)
    if do_sum:
        df = summarize_query_results(df, data, margin=margin)
    return(df)

def summarize_query_results(results, data, margin=True):
    results['Dx'] = data['Dx']
    results = results.groupby('Dx').sum().T
    if margin:
        categories = list(results.columns.categories) + ['All']
        ix = pd.CategoricalIndex(results.columns, categories=categories)
        results = results.reindex(columns=ix)
        results['All'] = results.sum(axis=1)
    return(results)
