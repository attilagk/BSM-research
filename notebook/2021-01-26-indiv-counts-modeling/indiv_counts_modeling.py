import pandas as pd
import numpy as np
from statsmodels.graphics import dotplots
from matplotlib import pyplot as plt

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

