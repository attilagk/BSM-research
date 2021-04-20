import pandas as pd
import numpy as np
import bsmutils

roadmap_rna_bname = bsmutils.get_bsmdir() + '/resources/roadmap-epigenomics/rna/expression/57epigenomes.'
proteinatlas_rna_bname = bsmutils.get_bsmdir() + '/resources/proteinatlas/expression/tissue_category_rna_brain_'

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

def gwascat_extract_genes(gwas, genecol='REPORTED GENE(S)', sep=', ', altsep=' - '):
    trait = gwas['DISEASE/TRAIT'].unique()
    s = gwas[genecol].dropna().str.replace(altsep, sep)
    s = set(s.str.split(sep).sum())
    return(s)

def gwascat_extract_reported_mapped_genes(gwas, trait):
    genecols = ['REPORTED GENE(S)', 'MAPPED_GENE']
    reported, mapped = (gwascat_extract_genes(gwas, genecol=genecol) for genecol in genecols)
    val = (reported, mapped)
    reported_mapped = reported.union(mapped)
    val = {trait + ' reported': reported, trait: reported_mapped}
    return(val)

def gwascat_multi_genesets(gwas, gwasPMID):
    D = {}
    def foo(trait='ADHD'):
            pmid = gwasPMID.loc[trait, 'PMID']
            d = gwascat_extract_reported_mapped_genes(gwas.loc[gwas['PUBMEDID'] == pmid], trait)
            D.update(d)
            return(d)
    for trait in gwasPMID.index:
        d = foo(trait)
    return(D)

