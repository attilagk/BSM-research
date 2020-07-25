import pandas as pd

def expandCNV(CNV, hgcb, ix):
    row = CNV.iloc[[ix], :]
    chromosome = row.loc[row.index[0], 'chromosome']
    start_cb = row.loc[row.index[0], 'start cytoband']
    end_cb = row.loc[row.index[0], 'end cytoband']
    hgcb = hgcb.loc[hgcb['chromosome'] == chromosome, :]
    name = list(hgcb['name'])
    val = range(name.index(end_cb), name.index(start_cb))
    return(val)
    hgcb = hgcb.iloc[range(name.index(end_cb), name.index(start_cb)), :]
    return(hgcb)
