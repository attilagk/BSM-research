import pandas as pd

def filtered_Dx_count(annot, filt='coding nonsyn', gb='Dx'):
    data = annot.copy()
    data[gb] = pd.Categorical(data[gb])
    data[gb] = data[gb].cat.set_categories(list(data[gb].cat.categories) + ['Total'])
    g = data.loc[data[filt].astype('bool'), [gb]].groupby(gb)
    s = g.size()
    s['Total'] = s.sum()
    df = pd.DataFrame(s, columns=[filt]).T
    return(df)

def filtered_Dx_counts(annot, filtl, gb='Dx'):
    dfl = [filtered_Dx_count(annot, y, gb='Dx') for y in filtl]
    df = pd.concat(dfl, axis=0)
    return(df)
