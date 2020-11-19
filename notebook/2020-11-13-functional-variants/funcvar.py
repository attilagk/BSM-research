import pandas as pd
import itertools
from bsmcalls import SNPnexus

clozukpath = '/home/attila/projects/bsm/resources/CLOZUK/supp-table-4.csv'

def filtered_Dx_count(annot, filt='coding nonsyn', gb='Dx', anygwas=None):
    data = annot.copy()
    data[gb] = pd.Categorical(data[gb])
    data[gb] = data[gb].cat.set_categories(list(data[gb].cat.categories) + ['Total'])
    filt_pass = data[filt].astype('bool')
    if anygwas is not None:
        filt_pass = filt_pass & anygwas
    g = data.loc[filt_pass, [gb]].groupby(gb)
    s = g.size()
    s['Total'] = s.sum()
    df = pd.DataFrame(s, columns=[filt]).T
    return(df)

def filtered_Dx_counts(annot, filtl, gb='Dx', anygwas=None):
    dfl = [filtered_Dx_count(annot, y, gb='Dx', anygwas=anygwas) for y in filtl]
    df = pd.concat(dfl, axis=0)
    return(df)

def get_geneset(df=pd.read_csv(clozukpath, skiprows=7), col='Gene(s) tagged'):
    ll = SNPnexus.str2set_setvalued(df, col, nonestr='', sepstr=',', listval=True)
    geneset = set(list(itertools.chain(*ll)))
    return(geneset)
