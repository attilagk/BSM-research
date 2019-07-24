import pandas as pd
import subprocess
import tempfile

def process_yifans_table(validation_path='/home/attila/projects/bsm/bsm-network/all_validations_190309_yifan.txt',
        outbn='/home/attila/projects/bsm/results/2019-07-23-commonsample-precrecall/somvar'):
    df = pd.read_csv(validation_path, delimiter='\t')
    df = df.loc[df['Manual Check'] == 'PASS', ['chrm', 'pos', 'ref', 'alt']]
    df.columns = ['#CHROM', 'POS', 'REF', 'ALT']
    df['#CHROM'] = df['#CHROM'].astype('category')
    df['#CHROM'] = df['#CHROM'].cat.set_categories([str(c) for c in range(22)] + ['X', 'Y'])
    df = df.sort_values(['#CHROM', 'POS'])
    somvar_path = outbn + '.tsv'
    df.to_csv(somvar_path, sep='\t', header=True, index=False)
    return(df)


def region_filter_callset(callset_fpath, somvar_fpath=None, somvar_df=None):
    regions = somvar_df[['#CHROM', 'POS']]
    regtmp = tempfile.NamedTemporaryFile()
    somvar_reg_fpath = regtmp.name
    regions.to_csv(somvar_reg_fpath, sep='\t', header=False, index=False)
    args0 = ['bcftools', 'query', '-R', somvar_reg_fpath, '-f', '%CHROM\t%POS\t%REF\t%ALT\n', callset_fpath]
    p0 = subprocess.Popen(args0, stdout=subprocess.PIPE)
    cnames = ['#CHROM', 'POS', 'REF', 'ALT']
    df = pd.read_csv(p0.stdout, sep='\t', header=None, names=cnames)
    regtmp.close()
    return(df)


def generate_record_id(vardf):
    row = vardf.iloc[0, :]
    ID = ':'.join([str(z) for z in row])
    ID = vardf.apply(lambda row: ':'.join([str(z) for z in row]), axis=1)
    return(ID)
