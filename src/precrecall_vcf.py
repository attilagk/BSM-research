import pandas as pd

def process_yifans_table(validation_path='/home/attila/projects/bsm/bsm-network/all_validations_190309_yifan.txt',
        outbn='/home/attila/projects/bsm/results/2019-07-23-commonsample-precrecall/somvar'):
    df = pd.read_csv(validation_path, delimiter='\t')
    df = df.loc[df['Manual Check'] == 'PASS', ['chrm', 'pos', 'ref', 'alt']]
    df.columns = ['#CHROM', 'POS', 'REF', 'ALT']
    df['#CHROM'] = df['#CHROM'].astype('category')
    df['#CHROM'] = df['#CHROM'].cat.set_categories([str(c) for c in range(22)] + ['X', 'Y'])
    df = df.sort_values(['#CHROM', 'POS'])
    regions = df[['#CHROM', 'POS']]
    somvar_path = outbn + '.tsv'
    df.to_csv(somvar_path, sep='\t', header=True, index=False)
    regions_path = outbn + '.regions'
    regions.to_csv(regions_path, sep='\t', header=False, index=False)
    return(df)
