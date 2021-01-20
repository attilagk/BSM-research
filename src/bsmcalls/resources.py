import pandas as pd
import numpy as np

roadmap_rna_bname = '/home/attila/projects/bsm/resources/roadmap-epigenomics/rna/expression/57epigenomes.'
proteinatlas_rna_bname = '/home/attila/projects/bsm/resources/proteinatlas/expression/tissue_category_rna_brain_'

def read_roadmap_rna(kind='RPKM', sampledict={'E071': 'BRN.HIPP.MID', 'E082': 'BRN.FET.F'}, suffix=False):
    if suffix:
        sampledict = dict(zip(sampledict.keys(), [x + '_' + kind for x in sampledict.values()]))
    fpath = roadmap_rna_bname + kind + '.pc'
    df = pd.read_csv(fpath, sep='\t', index_col=0, usecols=sampledict.keys())
    df = df.rename(sampledict, axis=1)
    return(df)

def read_roadmap_rna_RPKM_N(sampledict={'E071': 'BRN.HIPP.MID', 'E082': 'BRN.FET.F'}):
    l = [read_roadmap_rna(k, sampledict, suffix=True) for k in ['RPKM', 'N']]
    df = pd.concat(l, axis=1)
    return(df)

def read_proteinatlas_rna_brain(kind='elevated', index_col='Ensembl', usecols=['Gene', 'Gene synonym']):
    fpath = proteinatlas_rna_bname + kind + '.tsv'
    df = pd.read_csv(fpath, sep='\t', usecols=usecols + [index_col])
    df = df.set_index(index_col)
    if 'Gene synonym' in usecols:
        df['Gene all names'] = df['Gene synonym'].copy()
        df['Gene all names'] = df['Gene all names'].dropna().str.split(', ').apply(lambda x: set(x))
        df.dropna().apply(lambda x: x['Gene all names'].update({x['Gene']}), axis=1)
        df['Gene synonym'] = df['Gene synonym'].dropna().str.split(', ').apply(lambda x: set(x))
    return(df)

