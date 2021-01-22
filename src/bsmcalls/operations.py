import pandas as pd
import numpy as np
import functools
from bsmcalls import individuals
from scipy import stats

_level_names = ['Feature', 'Query']

def series_of_sets_intersect(series_of_sets, query):
    '''
    Intersect query set with a pandas.Series of sets
    '''
    if not isinstance(query, set):
        query = {query}
    s = series_of_sets.dropna()
    s = s.apply(lambda x: len(x.intersection(query))).astype(np.int16)
    s = s.reindex_like(series_of_sets).fillna(0, downcast=None).astype(np.int16)
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

features = {'ensembl_Gene': lambda x: x,
            'ensembl_Symbol': lambda x: dicts2sets(x, False),
            'ensembl_Predicted Function': lambda x: x,
            'near_gens_Overlapped Gene': lambda x: x,
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
    '''
    Query a feature in a dataframe

    Arguments
    queryitems: one or more strings to search for; see Examples 
    feature: a column name in data
    data: a DataFrame with rows: calls and columns: features

    Value: a DataFrame of one count per queryitem and call

    Details:
    Depending on the type and content of queryitems the output DataFrame has
    one or more columns with various semantics.  See Examples below

    Examples:

    I. Queryitems in a list

    > query(['coding nonsyn', 'missense'], feature='near_gens_Annotation', data=annot)

    Feature                                near_gens_Annotation
    Query                                         coding nonsyn missense
    Individual ID CHROM POS       Mutation
    CMC_MSSM_027  1     11973569  C/T                         0        0
                        67324428  A/T                         0        0
                        182008461 C/T                         0        0
                        207598957 C/T                         0        0
                        219342012 C/T                         0        0
    ...                                                     ...      ...

    II. Queryitems in a set valued dict

    > operations.query({'SCZ GWAS genes': gwasgenes}, 'near_gens_Overlapped Gene', data)

    Feature                                near_gens_Overlapped Gene
    Query                                             SCZ GWAS genes
    Individual ID CHROM POS       Mutation
    CMC_MSSM_027  1     11973569  C/T                              0
                        67324428  A/T                              0
                        182008461 C/T                              0
                        207598957 C/T                              0
                        219342012 C/T                              0
    ...                                                          ...

    III. When we search for any non NaN value in feature

    > operations.query(None, 'phast_Score', data)

    Feature                                phast_Score
    Query                                          any
    Individual ID CHROM POS       Mutation
    CMC_MSSM_027  1     11973569  C/T            False
                        67324428  A/T            False
                        182008461 C/T            False
                        207598957 C/T            False
                        219342012 C/T            False
    ...                                            ...
    UMB914        8     137927145 C/T            False
                  9     120507425 C/T            False
                        137710239 C/G            False
                  X     99255811  G/T            False
                        130133640 C/T             True

    [6426 rows x 1 columns]
    '''
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


def multiquery(querydict, data, do_sum=False, do_sort=False):
    '''
    Query multiple features at once

    Arguments
    querydict: a dict where keys are features and values are lists or dicts (see Example)
    data: an annot DataFrame
    do_sum: whether to summarize results by calling summarize_query_results
    do_sort: unsorted pandas multiindex gives a warning about low performance

    Value:
    a DataFrame of counts either for each call if do_sum=False or aggregated
    across all calls if do_sum=True

    Example:

    > querydict = {'sift_Prediction': ['Deleterious', 'Deleterious - Low Confidence'],
                   'tfbs_TFBS Name': None,
                   'near_gens_Overlapped Gene': {'SCZ GWAS genes': gwasgenes},
                   'ensembl_Gene': {'brain enriched HPA': set(pa_enriched.index), 'brain elevated HPA': set(pa_elevated.index)},
                   }
    > multiquery(querydict, data, do_sum=True, do_sort=False)

    Dx                                                      Control  SCZ  ASD  All
    Feature                   Query
    sift_Prediction           Deleterious                         8    9    6   23
                              Deleterious - Low Confidence        1    1    2    4
    tfbs_TFBS Name            any                                48   78   50  176
    near_gens_Overlapped Gene SCZ GWAS genes                     20   66   42  128
    ensembl_Gene              brain enriched HPA                 42   69   59  170
                              brain elevated HPA                211  344  280  835
    '''
    l = [query(v, k, data) for k, v in querydict.items()]
    df = pd.concat(l, axis=1)
    if do_sort:
        df = df.sort_index(axis=1)
    if do_sum:
        df = summarize_query_results(df, data)
    return(df)

def summarize_query_results(results, data, aggfun=None, chisq=True, margin=True):
    results['Dx'] = data['Dx']
    results = results.groupby('Dx')
    if aggfun is not None:
        results = results.apply(lambda x: x.groupby('Individual ID').sum().agg(aggfun)).T
    else:
        results = results.sum().T.astype('int64')
    if margin:
        categories = list(results.columns.categories) + ['All']
        ix = pd.CategoricalIndex(results.columns, categories=categories)
        results = results.reindex(columns=ix)
        results['All'] = results.sum(axis=1)
    if chisq:
        nsamples = individuals.get_nsamples(data, margin=True)
        results = chisquare_summary(results, nsamples)
    return(results)

def summarize_query_mean_sem(results, data):
    fundict = {'mean': np.mean, 'sem': lambda x: np.std(x) / (len(x) - 1)}
    df = summarize_query_results(results, data, aggfun=fundict.values(), chisq=False, margin=False)
    df = df.rename(columns={'<lambda>': 'sem'})
    return(df)

def chisquare_summary_row(observed, nsamples, onlyp=True):
    if 'All' not in observed.index:
        observed.loc['All'] = observed.sum()
    nsamples = pd.Series(nsamples) if isinstance(nsamples, dict) else nsamples
    expected = nsamples.to_numpy() / nsamples.loc['All'] * observed.loc['All']
    expected = pd.Series(expected, index=nsamples.index)
    val = stats.chisquare(observed.to_numpy()[:-1], expected.to_numpy()[:-1])
    if onlyp:
        val = val[1]
    return(val)

def chisquare_summary(summary, nsamples):
    summary = summary.copy()
    summary.columns = list(summary.columns)
    s = summary.apply(lambda x: chisquare_summary_row(x, nsamples, onlyp=False), axis=1)
    summary['chisq stat'] = s.apply(lambda x: x[0])
    summary['chisq p'] = s.apply(lambda x: x[1])
    return(summary)
