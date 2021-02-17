import pandas as pd
import numpy as np
from statsmodels.graphics import dotplots
from matplotlib import pyplot as plt
import statsmodels.api as sm
import patsy

def big_plot_matrix(responses, covariates, Dxcol, dropASD=False):
    covariates = covariates.select_dtypes(exclude='category')
    if dropASD:
        responses = responses.loc[Dxcol != 'C2']
        covariates = covariates.loc[Dxcol != 'C2']
        Dxcol = Dxcol[Dxcol != 'C2']
    nresp = responses.shape[1]
    ncovar = covariates.shape[1]
    fig, ax = plt.subplots(nresp, ncovar, sharey=False, figsize=(ncovar * 2, nresp * 2))
    for i, row in zip(range(nresp), responses.columns):
        response = responses[row]
        ax[i, 0].set_ylabel(row, rotation='horizontal', horizontalalignment='right')
        for j, col in zip(range(ncovar), covariates.columns):
            if i == 0:
                ax[i, j].set_title(col)
            if i == nresp - 1:
                ax[i, j].set_xlabel(col)
            ax[i, j].scatter(y=response, x=covariates[col], marker='|', c=Dxcol)
    return((fig, ax))

def my_dotplot(feature, mods):
    tvalues = pd.Series({endog: mods[endog].tvalues[feature] for endog in mods.keys()})
    pvalues = pd.Series({endog: mods[endog].pvalues[feature] for endog in mods.keys()})
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    fig.suptitle(feature + ' significance')
    ax[0].plot([0, 0], [0, 12])
    g = dotplots.dot_plot(tvalues, lines=tvalues.index, ax=ax[0])
    ax[0].set_xlim(np.array([-1, 1]) * 1.1 * tvalues.abs().max())
    ax[0].set_xlabel('t-value')
    g = dotplots.dot_plot(pvalues, lines=pvalues.index, ax=ax[1], show_names='right')
    ax[1].set_xscale('log')
    ax[1].set_xlabel('p-value')
    return((fig, ax))

def endog_binomial(feature, fitdata, proportion=False):
    success = fitdata[feature]
    if proportion:
        prop = success / fitdata['ncalls']
        return(prop)
    failure = fitdata['ncalls'] - success
    complement = 'NOT_' + feature
    df = pd.DataFrame({feature: success, complement: failure})
    return(df)

def my_logistic_fits(fitdata, endogname, exognames=['1', 'Dx', 'ageOfDeath', 'Dataset', 'AF', 'DP']):
    y = endog_binomial(endogname, fitdata, proportion=False)
    def helper(exogname):
        formula = ' + '.join(exognames[:exognames.index(exogname) + 1])
        X = patsy.dmatrix(formula, data=fitdata, return_type='dataframe')
        mod = sm.GLM(endog=y, exog=X, family=sm.families.Binomial()).fit()
        return((formula, mod))
    mods = dict([helper(exogname) for exogname in exognames])
    return(mods)
    formulas = [' + '.join(exognames[:exognames.index(x) + 1]) for x in exognames]
    Xs = [patsy.dmatrix(f, data=fitdata, return_type='dataframe') for f in formulas]
    mods = [sm.GLM(endog=y, exog=X, family=sm.families.Binomial()).fit() for X in Xs]
    return(Xs)
    X = patsy.dmatrix(formula, data=fitdata, return_type='dataframe')
    mod = sm.GLM(endog=y, exog=X, family=sm.families.Binomial()).fit()
