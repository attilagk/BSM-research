import pandas as pd
import os.path
import subprocess
import tempfile
import multi_prec_recall

def process_yifans_table(validation_path='/home/attila/projects/bsm/bsm-network/all_validations_190309_yifan.txt',
        outbn='/home/attila/projects/bsm/results/2019-07-23-commonsample-precrecall/somvar'):
    df = pd.read_csv(validation_path, delimiter='\t')
    df = df.loc[df['Manual Check'] == 'PASS', ['chrm', 'pos', 'ref', 'alt']]
    df.columns = ['#CHROM', 'POS', 'REF', 'ALT']
    df['#CHROM'] = df['#CHROM'].astype('category')
    df['#CHROM'] = df['#CHROM'].cat.set_categories([str(c) for c in range(22)] + ['X', 'Y'])
    df = df.sort_values(['#CHROM', 'POS'])
    var_tsvpath = outbn + '.tsv'
    df.to_csv(var_tsvpath, sep='\t', header=True, index=False)
    return(df)


def region_filter_callset(callset_vcfpath, var_tsvpath=None, var_df=None, PASS=True):
    '''
    Filters callset given the "regions" contained in var_df or var_tsvpath
    '''
    if var_df is None:
        var_df = pd.read_csv(var_tsvpath, delimiter='\t')
    regions = var_df[['#CHROM', 'POS']]
    regtmp = tempfile.NamedTemporaryFile()
    somvar_reg_fpath = regtmp.name
    regions.to_csv(somvar_reg_fpath, sep='\t', header=False, index=False)
    args0 = ['bcftools', 'query', '-R', somvar_reg_fpath, '-f',
            '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n', callset_vcfpath]
    p0 = subprocess.Popen(args0, stdout=subprocess.PIPE)
    cnames = ['#CHROM', 'POS', 'REF', 'ALT', 'FILTER']
    df = pd.read_csv(p0.stdout, sep='\t', header=None, names=cnames)
    regtmp.close()
    if PASS:
        df = df.loc[df['FILTER'] == 'PASS', :]
    return(df)


def generate_record_id(var_df):
    row = var_df.iloc[0, :]
    ID = ':'.join([str(z) for z in row])
    ID = var_df.apply(lambda row: ':'.join([str(z) for z in row]), axis=1)
    return(ID)


def prec_recall(callset_vcfpath,
        var_tsvpath='/big/results/bsm/2019-07-23-commonsample-precrecall/somvar.tsv',
        PASS=True):
    n_calls = multi_prec_recall.nrecords_in_vcf(callset_vcfpath, PASS=PASS)
    truecalls = region_filter_callset(callset_vcfpath=callset_vcfpath,
            var_tsvpath=var_tsvpath, PASS=PASS)
    caller = os.path.basename(callset_vcfpath)
    n_truecalls = len(truecalls)
    var_df = pd.read_csv(var_tsvpath, delimiter='\t')
    n_variants = len(var_df)
    #precision = n_truecalls / n_calls
    #recall = n_truecalls / n_variants
    #d = {'caller': caller, 'PASS': PASS, 'n_variants': n_variants, 'n_calls': n_calls, 'n_truecalls':n_truecalls, 'precision': precision, 'recall': recall}
    d = {'caller': caller, 'PASS': PASS, 'n_variants': n_variants, 'n_calls': n_calls, 'n_truecalls':n_truecalls}
    df = pd.DataFrame(d, index=[0])
    return(df)


#callers=['lofreqSomatic', 'somaticSniper', 'strelka2Germline', 'strelka2Somatic', 'Tnseq']
def prec_recall_all(callers=['lofreqSomatic', 'somaticSniper', 'strelka2Germline', 'strelka2Somatic', 'Tnseq'],
        callsetdir='/big/results/bsm/2018-02-22-ref-tissue-proj-testdata/wgs/vcf/snvs/',
        var_tsvpath='/big/results/bsm/2019-07-23-commonsample-precrecall/somvar.tsv'):
    callset_vcfpaths = [callsetdir + os.path.sep + c + '.vcf.gz' for c in callers]
    l_True = [prec_recall(callset_vcfpath=c, var_tsvpath=var_tsvpath,
        PASS=True) for c in callset_vcfpaths]
    l_False = [prec_recall(callset_vcfpath=c, var_tsvpath=var_tsvpath,
        PASS=False) for c in callset_vcfpaths]
    df = pd.concat(l_True + l_False)
    return(df)
