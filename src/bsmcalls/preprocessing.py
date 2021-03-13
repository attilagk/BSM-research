import pandas as pd
import numpy as np
import statsmodels as sm
import re


def convert2categorical(data, max_ncat=5):
    def helper(y):
        # if y contains only 1 kind of non-nan value
        if len(set(y.dropna())) < 2:
            return(y)
        # if y is a numeric categorical variable
        if not set(y.dropna()) - set(range(-1, max_ncat + 1)):
            max_level = np.int(max(set(y.dropna())))
            y = pd.Categorical(y, categories=range(max_level + 1), ordered=True)
        # if y is type object and contains at least one non-unique value
        elif y.dtype == 'object' and y[y == y.mode().values[0]].count() > 1:
            y = pd.Categorical(y)
        return(y)
    data = data.apply(helper, axis=0)
    if 'Dx' in data.columns:
        data['Dx'] = pd.Categorical(data['Dx'], categories=['Control', 'SCZ', 'ASD'], ordered=True)
    return(data)

def drop_unused_categories(data):
    val = data.copy()
    categorical = val.select_dtypes(include='category')
    for col in categorical.columns:
        val[col] = val[col].cat.remove_unused_categories()
    return(val)

def impute_vars(data, dropthrs=0.10):
    '''
    Impute variables with their mean or mode for numeric or categorical vars, respectively

    Arguments
    data: a data frame
    dropthrs: if the fraction of missing values exceeds it the variable is dropped

    Value: the imputed data frame
    '''
    # drop variables depending on dropthrs
    droplessthan = len(data) - len(data) // (1 / dropthrs)
    nobs = data.count()
    columns2drop = data.count().index[data.count() < droplessthan]
    data = data.drop(columns=columns2drop)
    # determine which columns to impute and how
    max_obs = len(data)
    cols2impute = set(nobs[nobs < max_obs].index)
    numeric_cols = set(data.select_dtypes(include=['float64', 'int64']).columns)
    categ_cols = set(data.select_dtypes(include=['category', 'object']).columns)
    # helper functions calculating the fillin value
    def impute_numeric(col): return(data[col].mean())
    def impute_categ(col): return(data[col].mode().values[0])
    # Create a copy so that the input 'data' isn't modified
    impdata = data.copy()
    # perform imputation separately for numeric and categorical columns
    for col in numeric_cols:
        impdata[col].fillna(value=impute_numeric(col), inplace=True)
    for col in categ_cols:
        impdata[col].fillna(value=impute_categ(col), inplace=True)
    return(impdata)

def prettify_colnames(data, repl='', pattern='[ ./\():	]+'):
    def foo(x):
        val = '_'.join(x) if isinstance(x, tuple) else x
        return(val)
    res = data.copy()
    res.columns = [foo(x) for x in res.columns]
    res = res.rename(lambda y: re.sub(pattern, repl, y), axis='columns')
    return(res)

def dummify_var(data, vname='Dx'):
    '''
    Dummify variable, typically Dx
    '''
    s = data[vname].astype('category')
    if len(s.cat.categories) != 2:
        raise TypeError(vname + ' is not a binary variable')
    data[vname] = s.cat.rename_categories([0, 1]).astype('int32')
    return(data)

def standardize_numvars(data):
    std_data = data.apply(lambda y: (y - y.mean()) / y.std()
                          if y.dtype == 'float64' or y.dtype == 'int64'
                          else y, axis=0)
    return(std_data)

def preprocess(data, impute=True, prettify=True, dummify=True, standardize=True):
    '''
    Preprocess data by imputing, prettifying, dummifying (Dx) and standardizing
    '''
    if impute:
        data = impute_vars(data)
    if prettify:
        data = prettify_colnames(data)
    if dummify:
        if 'Dx' in data.columns:
            data = dummify_var(data, vname='Dx')
    if standardize:
        data = standardize_numvars(data)
    return(data)

def dummify_df(df):
    numeric = df.select_dtypes(exclude='category')
    categorical = df.select_dtypes(include='category')
    dummy = pd.get_dummies(categorical)
    dfnum = pd.concat([numeric, dummy], axis=1)
    return(dfnum)

def collapse_categories(series, categdict):
    '''
    Collapse lists of categories into single categories

    Arguments
    series: a pandas Series of category dtype
    categdict: a dictionary whose keys are new categories and values are lists of old categories

    Value: the same series with the collapsed categories
    '''
    # take a single item of categdict
    newcategory, oldcategories = categdict.popitem()
    # remove old categories
    series = series.cat.remove_categories(oldcategories)
    # add new category and fill with it data points of old categories
    series = series.cat.add_categories(newcategory).fillna(newcategory)
    if len(categdict) == 0:
        return(series)
    else:
        return(collapse_categories(series, categdict))

def drop_category(data, categories='ASD', col='Dx'):
    s = data.copy()[col].cat.remove_categories(categories).dropna()
    data = data.copy().loc[s.index]
    data[col] = s
    return(data)
