import pandas as pd
import numpy as np
from scipy import stats

def cont_table(results, col, Dxcol=('Dx',)):
    Dx = results[Dxcol]
    Dx.name = 'Dx'
    feature = results[col]
    #feature.name = col[0]
    #if col[1] != 'any':
    #    feature.name += '_' + col[1]
    contingency = pd.crosstab(Dx, feature, margins=True, normalize=True)
    return(contingency)
