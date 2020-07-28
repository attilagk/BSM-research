import pandas as pd
import numpy as np

def expandCNVrow(CNV, hgcb, corr_cb, ix):
    row = CNV.iloc[[ix], :]
    chromosome = row.loc[row.index[0], 'chrom']
    start_cb = row.loc[row.index[0], 'start cytoband']
    # correct cytoband name if necessary to make it compatible with hgcb
    if (chromosome, start_cb) in corr_cb.index:
        start_cb = corr_cb.loc[[(chromosome, start_cb)], 'corr start'][0][1]
    end_cb = row.loc[row.index[0], 'end cytoband']
    if (chromosome, end_cb) in corr_cb.index:
        end_cb = corr_cb.loc[[(chromosome, end_cb)], 'corr end'][0][1]
    # get all cytobands for the given chromosome from hgcb
    hgcb = hgcb.loc[hgcb['chrom'] == chromosome, :]
    name = list(hgcb['name'])
    # find the indices of the start and end cytobands
    hgcb_ix = list((name.index(start_cb), name.index(end_cb)))
    hgcb_ix.sort()
    # Expand row by the number of cytobands defined by the start and end band.
    # When the start and end band are the same the expansion factor is trivially 1.
    row = row.iloc[np.zeros(np.diff(hgcb_ix) + 1, dtype=int), :]
    hgcb_chunk = hgcb.iloc[range(hgcb_ix[0], hgcb_ix[1] + 1), :]
    row.index = hgcb_chunk.index
    # Now take info from hgcb to the CNV row
    row['hg cytoband'] = hgcb_chunk['name']
    row['hg chr cytoband'] = [(c, n) for c, n in zip(hgcb_chunk['chrom'], hgcb_chunk['name'])]
    row['start'] = hgcb_chunk['start']
    row['end'] = hgcb_chunk['end']
    row.drop(columns=['start cytoband', 'end cytoband'], inplace=True)
    return(row)


def expandCNV(CNV, hgcb, corr_cb):
    '''
    Expands the table of schizophrenia related CNVs from possibly multiple
    cytogenic band to single band rows

    Arguments
    CNV: info on CNVs preprocessed in 2020-07-24-szdb.ipynb
    hgcb: mapping between cytogenic bands to base positions
    corr_cb: mapping between correct and incorrect cytogenic band names manually created in 2020-07-24-szdb.ipynb

    Value: the expanded CNV
    '''
    l = [expandCNVrow(CNV, hgcb, corr_cb, ix) for ix in range(len(CNV))]
    CNV = pd.concat(l, axis=0)
    return(CNV)


def bedifyCNV(CNV, hgcb):
    '''
    Count occurrence of each cytogenic band in szdb expanded CNV table and then bedify

    Arguments
    CNV: expanded CNV returned by expandCNV
    hgcb: as above

    Value: the bedified CNV
    '''
    gb = CNV.groupby(by='hg chr cytoband')
    first = gb.first()[['start', 'end']]
    chrom = pd.DataFrame({'chrom': [y[0] for y in first.index]})
    chrom.index = first.index
    count = pd.DataFrame({'count': gb.size()})
    present = pd.concat([chrom, first, count], axis=1)
    # some odds and ends for present
    present.index = pd.MultiIndex.from_tuples(present.index)
    # add missing zero count cyto bands
    missing = hgcb.loc[list(set(hgcb.index) - set(present.index)), ['chrom', 'start', 'end']]
    missing['count'] = 0
    val = pd.concat([present, missing])
    categories = [str(y) for y in range(1, 23)] + list('XY')
    val['chrom'] = pd.Categorical(val['chrom'], ordered=True, categories=categories)
    val.sort_values(by=['chrom', 'start'], inplace=True)
    return(val)
