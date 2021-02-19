import pandas as pd
import numpy as np
from statsmodels.graphics import dotplots
from matplotlib import pyplot as plt
import statsmodels.api as sm
import patsy
import scipy.stats as stats


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

def my_fits(fitdata, endogname, exognames=['1', 'Dx', 'ageOfDeath', 'Dataset', 'AF', 'DP'], family=sm.families.Poisson()):
    if isinstance(family, sm.families.Binomial):
        y = endog_binomial(endogname, fitdata, proportion=False)
    else:
        y = fitdata[endogname]
    def helper(exogname):
        formula = ' + '.join(exognames[:exognames.index(exogname) + 1])
        X = patsy.dmatrix(formula, data=fitdata, return_type='dataframe')
        mod = sm.GLM(endog=y, exog=X, family=family).fit()
        return((formula, mod))
    mods = dict([helper(exogname) for exogname in exognames])
    return(mods)

def r_star_residuals(mod):
    r_D = mod.resid_deviance
    r_P = mod.resid_pearson
    r_star = r_D + np.log(r_P / r_D) / r_D
    return(r_star)

def modsel_dotplot(mods):
    AIC = [mods[f].aic for f in mods.keys()]
    BIC = [mods[f].bic_llf for f in mods.keys()]
    def get_pvalues(m, param='Dx[T.SCZ]'):
        try:
            return(m.pvalues[param])
        except KeyError:
            return(np.nan)

    pvalue_SCZ = [get_pvalues(mods[f], param='Dx[T.SCZ]') for f in mods.keys()]
    pvalue_ASD = [get_pvalues(mods[f], param='Dx[T.ASD]') for f in mods.keys()]
    fig, ax = plt.subplots(1, 4, figsize=(12, 4))
    g = dotplots.dot_plot(AIC, lines=list(mods.keys()), ax=ax[0])
    g = dotplots.dot_plot(BIC, lines=list(mods.keys()), ax=ax[1], show_names='right')
    g = dotplots.dot_plot(pvalue_SCZ, lines=list(mods.keys()), ax=ax[2], show_names='right')
    g = dotplots.dot_plot(pvalue_ASD, lines=list(mods.keys()), ax=ax[3], show_names='right')
    ax[2].set_xscale('log')
    ax[3].set_xscale('log')
    ax[0].set_title('Model fit: AIC')
    ax[1].set_title('Model fit: BIC')
    ax[2].set_title('Dx[T.SCZ]: p-value')
    ax[3].set_title('Dx[T.ASD]: p-value')
    return((fig, ax))

def QQ_four_residual_types(mod):
    fig, ax = plt.subplots(2, 2, figsize=(12, 12))
    g = sm.qqplot(mod.resid_anscombe_scaled, stats.norm, line='45', ax=ax[0, 0])
    ax[0, 0].set_title('Scaled Anscombe residuals')
    g = sm.qqplot(mod.resid_deviance, stats.norm, line='45', ax=ax[0, 1])
    ax[0, 1].set_title('Deviance residuals')
    g = sm.qqplot(mod.resid_pearson, stats.norm, line='45', ax=ax[1, 0])
    ax[1, 0].set_title('Pearson residuals')
    r_star = r_star_residuals(mod)
    g = sm.qqplot(r_star, stats.norm, line='45', ax=ax[1, 1])
    ax[1, 1].set_title('$r^\star$ residuals')
    return((fig, ax))

def QQ_rstar_residual(models):
    fig, axi = plt.subplots(4, 3, sharex=True, sharey=True, figsize = (12, 16))
    for f, ax in zip(models.keys(), np.ravel(axi)):
        m = models[f]
        r_star = r_star_residuals(m)
        g = sm.qqplot(r_star, stats.norm, line='45', ax=ax)
        ax.set_title(f)
    return((fig, ax))
