import pandas as pd
import subprocess
import io
import bsmutils
from bsmcalls import operations

hs37d5_bed = '/big/data/refgenome/GRCh37/dna/hs37d5.bed'

def read_scz_gwas():
    fpath = bsmutils.get_bsmdir() + '/resources/CLOZUK/supp-table-4.csv'
    gwas = pd.read_csv(fpath, skiprows=7)
    gwas['Chromosome'] = pd.Categorical(gwas['Chromosome'], ordered=True, categories=[str(y) for y in range(1, 23)] + list('XY'))
    gwas.sort_values(['Chromosome', 'Start (BP)'], inplace=True)
    return(gwas)


def order_coordinates(df, coord_cols=['Chromosome', 'Start (BP)', 'End (BP)']):
    '''
    Order df according to genomic coordinates: chromosomes (1, 2, 3, ..., 22, X, Y) then start base
    '''
    categories = [str(y) for y in range(1, 23)] + list('XY')
    df[coord_cols[0]] = pd.Categorical(df[coord_cols[0]], ordered=True, categories=categories)
    df.sort_values(coord_cols[0:2], inplace=True)
    return(df)


def complement_intervals(df, coord_cols=['Chromosome', 'Start (BP)', 'End (BP)'], onebased=True):
    '''
    Complement genomic intervals in df to the whole genome

    Arguments
    df: the data frame; each row comes with a unique genomic interval
    coord_cols: names of columns giving genomic coordinates
    onebased (boolean): True if coordinates are one based; False if zero based

    Value: the data frame with the complement intervals added as new rows.

    Details:  The returned data frame has empty rows at the complement
    intervals.  In those rows coord_cols still have nonempty values, of
    course.
    '''
    df = order_coordinates(df, coord_cols)
    if onebased:
        df[coord_cols[1]] = df.loc[:, coord_cols[1]] - 1
    df_coord = df.loc[:, coord_cols]
    Bbedpath = '/tmp/order_coordinates_B.bed'
    df_coord.to_csv(Bbedpath, sep='\t', header=False, index=False)
    Abedpath = hs37d5_bed
    l = ['bedtools', 'subtract', '-a', Abedpath, '-b', Bbedpath]
    p = subprocess.run(l, capture_output=True)
    cdf = pd.read_csv(io.BytesIO(p.stdout), sep='\t', names=coord_cols) # cdf for complement df
    cdf = order_coordinates(cdf, coord_cols)
    val = pd.concat([df, cdf])
    val = order_coordinates(val, coord_cols)
    return(val)

def annotate_with_gwas_loci(data, gwas=read_scz_gwas(), modify_return_data=True):
    CHROM = data.index.get_level_values('CHROM')
    POS = data.index.get_level_values('POS')
    df = pd.DataFrame({'CHROM': CHROM, 'POS': POS}, index=data.index)
    selcols = ['Chromosome', 'Start (BP)', 'End (BP)', 'P-value', 'Length (KB)']
    for locus in gwas['Locus']:
        chrom, start, end, pval, length = gwas.loc[gwas['Locus'] == locus, selcols].to_numpy()[0]
        b = (df['CHROM'] == chrom) & (start <= df['POS']) & (df['POS'] <= end)
        df.loc[b, 'SCZ GWAS locus'] = locus 
        df.loc[b, 'SCZ GWAS p-value'] = pval 
        df.loc[b, 'SCZ GWAS length (KB)'] = length 
    df['SCZ GWAS locus'] = df['SCZ GWAS locus'].astype('Int64')#.astype('str')
    coding = ~ data['near_gens_Overlapped Gene'].isna()
    s = df.loc[coding, 'SCZ GWAS locus'].dropna()
    df['SCZ GWAS locus, coding'] = s.reindex(data.index)
    df.loc[coding, 'SCZ GWAS p-value, coding'] = df.loc[coding, 'SCZ GWAS p-value']
    if modify_return_data:
        newdata = pd.concat([data, df], axis=1)
        return(newdata)
    return(df)
