import pandas as pd
import subprocess
import io

hs37d5_bed = '/big/data/refgenome/GRCh37/dna/hs37d5.bed'

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
