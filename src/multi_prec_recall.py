import pandas as pd
import numpy as np
import os
import os.path
import tempfile
import shutil
import subprocess
import glob
import truth_sets_aaf as tsa
import seaborn
import matplotlib.pyplot as plt

__callsets__ = ['strelka2Germline', 'strelka2Germline2s', 'MosaicForecast',
        'lofreqSomatic', 'MuTect2', 'strelka2Somatic', 'somaticSniper']
__callsetmaindir__ = '/home/attila/projects/bsm/results/calls/mixing-experiment/'
__truthsetmaindir__ = '/home/attila/projects/bsm/results/2019-03-18-truth-sets/'
__outmaindir__ = '/home/attila/projects/bsm/results/2019-08-15-benchmark-calls/'
__expmsubdir__ = 'truthset/aaf/exp_model/lambda_'
__expmsubdir1__ = 'filtered4exp_model/lambda_'
__addthreads__ = '7'
__allthreads__ = str(int(__addthreads__) + 1)
__markers__ = ['o', 'X', 's', 'P', 'd', '^', 'v']


def getVCFpaths(callsetbn=None, region='chr22', vartype='snp', lam='0.04',
        s2g='-2', case_sample='mix1', control_sample='mix3', callsetdir=None):
    '''
    Create pathname for various input, output, intermediate VCF files.

    Parameters:
    callsetbn: None, 'Tnseq.vcf.gz' or ['lofreqSomatic.vcf.gz', 'Tnseq.vcf.gz']
    region: chr22, chr1_2 or autosomes
    vartype: snp or indel
    lam: '0.04' or '0.2'
    s2g: '-2', '-3' or '-4'
    case_sample: 'mix1', 'mix2' or 'mix3'
    control_sample: 'mix1', 'mix2', 'mix3' or 'no_ctr'

    Returns:
    
    a dictionary of pathnames

    The keys of the dictionary are as follows:

    callset: the original callset without filtering

    prepared_callset_dir: directory of VCFs filtered for region and vartype
    
    prepared_callset: filtered for region and vartype
    
    reduced_truthset: the truthset according to the exp_model with parameters
    lam, s2g, case_sample
    
    discarded_from_truthset: the complementer set of reduced_truthset relative to
    the original truthset based on the original mixes

    reduced_callset: created from reduced_truthset by removing the
    "nonvariants" of discarded_from_truthset

    Details:

    If callsetbn is None then a list of all pathnames are used that match the
    pattern '*.vcf.gz'.  If callsetbn is a list of file basenames then those
    will be extended into pathnames.  If callsetbn is a string of a single
    basename then that will be extended into a single pathname.
    '''
    # some pieces of pathnames
    sample_pair = case_sample + '-' + control_sample
    subdir0 = region + os.path.sep + vartype + os.path.sep
    subdir1 = sample_pair + os.path.sep + region + os.path.sep + vartype + os.path.sep
    subdir2 = lam + '/s2g_' + s2g + os.path.sep
    subdir3 = subdir2 + case_sample + os.path.sep
    if vartype == 'snp':
        alt_vartype = 'snvs'
    elif vartype == 'indel':
        alt_vartype = 'indels'
    # directories
    truthsetdir = __truthsetmaindir__ + subdir0 + __expmsubdir__ + subdir3
    if callsetdir is None:
        callsetdir = __callsetmaindir__ + sample_pair + os.path.sep + alt_vartype + os.path.sep
    prepared_callset_dir = __outmaindir__ + subdir1
    reduced_callset_dir = __outmaindir__ + subdir1 + __expmsubdir1__ + subdir2
    # filepaths
    reduced_truthset =  truthsetdir + 'complete.vcf.gz'
    discarded_from_truthset =  truthsetdir + 'discarded-complete.vcf.gz'
    if callsetbn is None:
        callset = glob.glob(callsetdir + '*.vcf.gz')
        prepared_callset = [prepared_callset_dir + os.path.basename(y) for y in callset]
        reduced_callset = [reduced_callset_dir + os.path.basename(y) for y in callset]
    elif isinstance(callsetbn, list):
        callset = [callsetdir + y for y in callsetbn]
        prepared_callset = [prepared_callset_dir + y for y in callsetbn]
        reduced_callset = [reduced_callset_dir + y for y in callsetbn]
    else:
        callset = callsetdir + callsetbn
        prepared_callset = prepared_callset_dir + callsetbn
        reduced_callset = reduced_callset_dir + callsetbn
    VCFpaths = {'reduced_truthset': reduced_truthset,
            'discarded_from_truthset': discarded_from_truthset, 'callset':
            callset, 'prepared_callset_dir': prepared_callset_dir,
            'prepared_callset': prepared_callset, 'reduced_callset_dir':
            reduced_callset_dir, 'reduced_callset': reduced_callset}
    return(VCFpaths)


def count_vcf_records(VCFpath):
    '''
    Counts records in VCF

    Parameter:
    VCFpath: the path to the VCF

    Returns: the count (integer)
    '''
    args1 = ['bcftools', 'view', '--threads', __addthreads__, '-H', VCFpath]
    args2 = ['wc', '-l']
    proc1 = subprocess.Popen(args1, shell=False, stdout=subprocess.PIPE)
    proc2 = subprocess.Popen(args2, shell=False, stdout=subprocess.PIPE, stdin=proc1.stdout)
    proc1.stdout.close()
    nrec = proc2.communicate()[0]
    nrec = int(nrec) # turn bytesliteral (e.g. b'5226\n') to integer
    return(nrec)


def get_callsetbn(vartype='snp', case_sample='mix1', control_sample='mix2', from_prepared_callset_dir=False):
    '''
    Gets callsetbn (see getVCFpaths) for available callsets for vartype, case_sample and control_sample

    Parameters:
    vartype: snp or indel
    case_sample: 'mix1', 'mix2' or 'mix3'
    control_sample: 'mix1', 'mix2' or 'mix3'

    Returns: the callsetbn (list)
    '''
    VCFpaths = getVCFpaths(vartype=vartype, case_sample=case_sample, control_sample=control_sample)
    if not from_prepared_callset_dir:
        callsets = VCFpaths['callset']
        callsetbn = [os.path.basename(y) for y in callsets if count_vcf_records(y) > 0]
    else:
        callsets = glob.glob(VCFpaths['prepared_callset_dir'] + '*.vcf.gz')
        callsetbn = [os.path.basename(y) for y in callsets]
    return(callsetbn)


def prepare4prec_recall(region='chr22', vartype='snp', case_sample='mix1',
        control_sample='mix2'):
    '''
    Run prepare4prec-recall shell script on initial callset VCFs for a given
    region and variant type

    Parameter(s):
    region: chr22, chr1_2 or autosomes
    vartype: snp or indel
    case_sample: 'mix1', 'mix2' or 'mix3'
    control_sample: 'mix1', 'mix2' or 'mix3'

    Returns: a list of the pathname of output VCFs
    '''
    callsetbn = get_callsetbn(vartype, case_sample=case_sample,
            control_sample=control_sample, from_prepared_callset_dir=False)
    VCFpaths = getVCFpaths(callsetbn=callsetbn, region=region,
            vartype=vartype, case_sample=case_sample, control_sample=control_sample)
    outdir = VCFpaths['prepared_callset_dir']
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    csetVCFs = VCFpaths['callset']
    # wrapper to the helper function
    def preparer(y):
        val = do_prepare4prec_recall(csetVCF=y, outdir=outdir, region=region,
                vartype=vartype, normalize=True, PASS=True)
        return(val)
    val = [preparer(y) for y in csetVCFs]
    return(val)


def do_prepare4prec_recall(csetVCF, outdir, region, vartype='snp',
        normalize=True, PASS=True, overwrite=False):
    '''
    Run prepare4prec-recall shell script on initial callset VCF for a given
    region and variant type with or without normalization and PASS filtering

    Parameter(s):
    csetVCF: path to the input callset
    region: chr22 or autosomes
    vartype: snp or indel
    normalize: (True|False) whether to perform normalizaton
    PASS: (True|False) whether to PASS filter

    Returns: the pathname of output VCF
    '''
    if normalize:
        normalize_opt = []
    else:
        normalize_opt = ['-n']
    if PASS:
        PASS_opt = ['-P']
    else:
        PASS_opt = []
    if region == 'autosomes':
        r_opt = []
    elif region == 'chr1_2':
        r_opt = ['-r1,2']
    elif region == 'chr22':
        r_opt = ['-r22']
    args1 = ['prepare4prec-recall', '-p', __allthreads__, '-v',
            vartype, '-Oz'] + normalize_opt + PASS_opt + r_opt + [csetVCF]
    outVCF = outdir + os.path.basename(csetVCF)
    args2 = ['bcftools', 'view', '--threads', __addthreads__, '-Oz', '-o', outVCF]
    if os.path.isfile(outVCF) and not overwrite:
        return(outVCF)
    proc1 = subprocess.Popen(args1, shell=False, stdout=subprocess.PIPE)
    proc2 = subprocess.run(args2, shell=False, stdout=subprocess.PIPE, stdin=proc1.stdout)
    args3 = ['bcftools', 'index', '--threads', __addthreads__, '--tbi', outVCF]
    proc3 = subprocess.run(args3)
    return(outVCF)


def reduce_prepared_callsets(callsetbn=None, region='chr22', vartype='snp', lam='0.04',
        s2g='-2', case_sample='mix1', control_sample='mix3', overwrite=False):
    '''
    Reduces (discards nonvariants from) the prepared callsets for a given region, vartype and exp_model

    Parameters:
    callsetbn: None or a list of callset VCF basenames like ['Tnseq.vcf.gz',...]
    region: chr22 or autosomes
    vartype: snp or indel
    lam: '0.04' or '0.2'
    s2g: '-2', '-3' or '-4'
    case_sample: 'mix1', 'mix2' or 'mix3'
    control_sample: 'mix1', 'mix2' or 'mix3'
    overwrite: whether to overwrite existing reduced callsets

    Returns: a list of pathnames of the reduced callsets
    '''
    VCFpaths = getVCFpaths(region=region, vartype=vartype, case_sample=case_sample, control_sample=control_sample)
    if callsetbn is None:
        callsetbn = get_callsetbn(vartype, case_sample=case_sample,
                control_sample=control_sample, from_prepared_callset_dir=False)
    VCFpaths = getVCFpaths(callsetbn=callsetbn, region=region,
            vartype=vartype, lam=lam, s2g=s2g, case_sample=case_sample, control_sample=control_sample)
    def helper(prepared_cset):
        discarded_tset = VCFpaths['discarded_from_truthset']
        red_cset_dir = VCFpaths['reduced_callset_dir']
        reduced_callset = red_cset_dir + os.path.basename(prepared_cset)
        if not overwrite and os.path.isfile(reduced_callset):
            return(reduced_callset)
        if not os.path.isdir(red_cset_dir):
            os.makedirs(red_cset_dir)
        reduced_callset_tbi = reduced_callset + '.tbi'
        tempdir = tempfile.TemporaryDirectory()
        tempVCF = tempdir.name + os.path.sep + '0000.vcf.gz'
        tempVCFtbi = tempVCF + '.tbi'
        args1 = ['bcftools', 'index', '--threads', __addthreads__, '--tbi', discarded_tset]
        subprocess.run(args1)
        args2 = ['bcftools', 'isec', '-C', '-Oz', '-p', tempdir.name,
                prepared_cset, discarded_tset]
        subprocess.run(args2)
        shutil.move(tempVCF, reduced_callset)
        shutil.move(tempVCFtbi, reduced_callset_tbi)
        return(reduced_callset)
    red_callsets = [helper(y) for y in VCFpaths['prepared_callset']]
    return(VCFpaths)


def prec_recall_one_truthset(truthset, callsets):
    args = ['prec-recall-vcf', '-p', __allthreads__, '-t', truthset] + callsets
    proc1 = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    prcsv = pd.read_csv(proc1.stdout)
    return(prcsv)


def reduce_precrecall(region='chr22', vartype='snp', lam='0.04',
        s2g='-2', case_sample='mix1', control_sample='mix2'):
    VCFpaths = reduce_prepared_callsets(region=region, vartype=vartype,
            lam=lam, s2g=s2g, case_sample=case_sample, control_sample=control_sample, overwrite=False)
    truthset = VCFpaths['reduced_truthset']
    callsets = VCFpaths['reduced_callset']
    pr = prec_recall_one_truthset(truthset=truthset, callsets=callsets)
    pr['region'] = region
    pr['vartype'] = vartype
    pr['lam'] = lam
    pr['s2g'] = s2g
    pr['case_sample'] = case_sample
    pr['control_sample'] = control_sample
    pr = pr_astype(pr, alphabetical=True)
    return(pr)


def vcf_exists(vartype, control_sample):
    VCFpaths = getVCFpaths(vartype=vartype, control_sample=control_sample)
    return(bool(len(VCFpaths['callset'])))


def prepare_reduce_precrecall(region='chr22', vartype='snp', case_sample='mix1'):
    '''
    Prepare and reduce callset and calculate precision and recall for a given
    region and variant type
    '''
    def process1exp_model(lam, s2g, control_sample):
        val = prepare4prec_recall(region=region, vartype=vartype,
                case_sample=case_sample, control_sample=control_sample)
        pr = reduce_precrecall(region=region, vartype=vartype, lam=lam,
                s2g=s2g, case_sample=case_sample, control_sample=control_sample)
        return(pr)
    lams = ['0.04', '0.2']
    s2gs = ['-2', '-3', '-4']
    control_samples = ['mix1', 'mix2', 'mix3', 'no_ctr']
    l = [process1exp_model(lam=l, s2g=g, control_sample=s) for l in lams for g in
            s2gs for s in control_samples if vcf_exists(vartype, s)]
    if len(l):
        pr = pd.concat(l)
        pr = pr_astype(pr, alphabetical=True)
        return(pr)
    else:
        return(None)


def vmc_prepare_reduce_precrecall(csetVCF, region='chr22', vartype='snp',
        machine='Ada'):
    '''
    Performs all three major steps: prepares the VCF, reduces it, and uses it
    for precision--recall calculations.
    '''
    callsetbn = [os.path.basename(csetVCF)]
    VCFpaths = getVCFpaths(callsetbn=callsetbn, region=region, vartype=vartype)
    outdir = VCFpaths['prepared_callset_dir']
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outVCF = do_prepare4prec_recall(csetVCF=csetVCF, outdir=outdir,
            region=region, vartype=vartype, normalize=False, PASS=False, overwrite=False)
    # model specific operations
    case_sample='mix1'
    control_sample='mix3'
    def helper(lam, s2g):
        VCFpaths = reduce_prepared_callsets(callsetbn=callsetbn, region=region, vartype=vartype, lam=lam,
                s2g=s2g, case_sample=case_sample, control_sample=control_sample, overwrite=False)
        csetVCF = VCFpaths['reduced_callset'][0]
        tsetVCF = VCFpaths['reduced_truthset']
        pr = vmc_precrecall(csetVCF=csetVCF, tsetVCF=tsetVCF)
        pr['region'] = region
        pr['vartype'] = vartype
        pr['lam'] = lam
        pr['s2g'] = s2g
        pr['case_sample'] = case_sample
        pr['control_sample'] = control_sample
        pr['machine'] = machine
        return(pr)
    lams = ['0.04', '0.2']
    s2gs = ['-2', '-3', '-4']
    l = [helper(lam=l, s2g=g) for l in lams for g in s2gs]
    pr = pd.concat(l)
    pr = pr_astype(pr, vmc_pr=True)
    return(pr)


#def run_all():
def prepare_reduce_precrecall_all():
    '''
    Prepare and reduce callset and calculate precision and recall for all
    regions and variant types
    '''
    regions = ['chr22', 'chr1_2', 'autosomes']
    vartypes = ['snp', 'indel']
    l = [prepare_reduce_precrecall(region=r, vartype=v) for r in regions for v
            in vartypes]
    l = [y for y in l if y is not None]
    pr = pd.concat(l)
    pr = pr_astype(pr, alphabetical=True)
    return(pr)


def read_runtime(fpath, region, machine):
    '''
    Read POSIX formatted runtimes into a pandas DataFrame
    '''
    rt = pd.read_csv(fpath, delim_whitespace=True, names=['type', 'runtime'])
    rt['region'] = region
    rt['machine'] = machine
    if region == 'chr1_2':
        length = 249250621 + 243199373
    if region == 'chr22':
        length = 51304566
    rt['region_length'] = length
    return(rt)


def correct_vmc_pr(pr, corr_f = (249250621 + 243199373) / 86682278 ):
    corr_pr = pr.copy()
    corr_pr['recall'] = pr['recall'] * corr_f
    return(corr_pr)


def pr_astype(pr, vmc_pr=False, alphabetical=False):
    '''
    Set data types for a precision recall data frame

    Parameters:
    pr: a precision recall data frame
    vmc_pr: True if pr is from VariantMetaCaller
    alphabetical: whether to sort category names alphabetically or according to the __callsets__ variable

    Returns: the data frame with the same data but corrected data types
    '''
    keys = ['region', 'vartype', 'lam', 's2g', 'case_sample', 'control_sample']
    if vmc_pr:
        pr['s2g'] = np.int64(pr['s2g']) # crucial for consistency
        keys = keys + ['machine', 'chrom', 'ref', 'alt']
    else:
        keys.append('callset')
    d = {k: 'category' for k in keys}
    val = pr.astype(d)
    if not vmc_pr:
        if alphabetical:
            categs = sorted(val['callset'].cat.categories, key=str.lower)
        else:
            categs = __callsets__
        val['callset'] = val['callset'].cat.set_categories(categs)
    return(val)


def read_pr_csv(csvpath, vmc_pr=False):
    '''
    Read precision recall data from a CSV into a data frame and set data types
    '''
    pr = pd.read_csv(csvpath)
    if not vmc_pr:
        pr = fix_names(pr)
    pr = pr_astype(pr, vmc_pr=vmc_pr)
    return(pr)


def replace_categ(df, column='callset', old='Tnseq', new='MuTect2'):
    df = df.copy()
    l = list(df['callset'])
    l = [x.replace(old, new) for x in l]
    df['callset'] = l
    return(df)

def replace_colname(df, old='s2g', new='s2g'):
    df = df.copy()
    l = [x.replace(old, new) for x in df.columns]
    df.columns = l
    return(df)


def fix_names(df):
    df = replace_categ(df, column='callset', old='Tnseq', new='MuTect2')
    df = replace_categ(df, column='callset', old='TNseq', new='MuTect2')
    return(df)


def singles2paireds(pr):
    '''
    Take rows for which control_sample is missing ('no_ctr') and add them to
    all rows with control sample.
    '''
    sel_rows = (pr['control_sample'] == 'no_ctr')
    def helper(control_sample='mix1'):
        df = pr.loc[sel_rows, :].copy()
        df['control_sample'] = control_sample
        return(df)
    csamples = ['mix1', 'mix2', 'mix3']
    l = [helper(control_sample=cs) for cs in csamples]
    res = pd.concat([pr] + l)
    return(res)


def plotter_vmc1(pr, vmc_pr, lam=0.2, region='chr1_2', s2g=-3,
        case_sample='mix1', control_sample='mix3', vartype='snp',
        callset=['strelka2Germline2s', 'strelka2Somatic', 'MuTect2', 'lofreqSomatic', 'somaticSniper']):
    '''
    Parameters:

    pr: a precision recall data frame
    vmc_pr: a precision recall data frame for VariantMetaCaller
    lam: '0.04' or '0.2'
    region: chr22 or autosomes
    s2g: '-2', '-3' or '-4'
    case_sample: 'mix1', 'mix2' or 'mix3'
    control_sample: 'mix1', 'mix2' or 'mix3'
    vartype: snp or indel
    callset: a list of caller names that were input to VMC

    Returns:
    a FacetGrid plot object; Precision-recall curve for VariantMetaCaller
    '''
    seaborn.set()
    seaborn.set_context('paper')
    sel_rows = (pr['case_sample'] == case_sample) & (pr['control_sample'] == control_sample) & \
                    (pr['s2g'] == s2g) & (pr['region'] == region) & (pr['lam'] == lam) & \
                                    (pr['callset'].isin(callset)) & (pr['vartype'] == vartype)
    pr_sset = pr.loc[sel_rows, :]
    machine = 'Ada'
    vmc_sel_rows = (vmc_pr['case_sample'] == case_sample) & (vmc_pr['control_sample'] == control_sample) & \
                    (vmc_pr['s2g'] == s2g) & (vmc_pr['region'] == region) & (vmc_pr['lam'] == lam) & \
                            (vmc_pr['machine'] == machine) & (vmc_pr['vartype'] == vartype)
    vmc_pr_sset = vmc_pr.loc[vmc_sel_rows, :]
    fg = seaborn.FacetGrid(data=pr_sset, margin_titles=True, aspect=1,
            hue='callset', sharey=True, hue_kws=dict(marker=__markers__))
    fg.axes[0][0].plot(vmc_pr_sset['recall'], vmc_pr_sset['precision'], color='black', linestyle='-')
    fg.axes[0][0].plot(vmc_pr_sset['recall'], vmc_pr_sset['precision_estim'], color='black', linestyle=':')
    fg = fg.map(plt.plot, 'recall', 'precision')
    fg = fg.add_legend()
    return(fg)


def plotter1b(pr, vmc_pr=None, sample='mix1', s2g=-2, vartype='snp'):
    '''
    Precision-recall plot; rows by s2g and columns by lambda
    '''
    seaborn.set()
    seaborn.set_context('talk')
    sel_rows = (pr['sample'] == sample) & (pr['s2g'] == s2g) & (pr['vartype'] == vartype)
    df_sset = fix_names(pr.loc[sel_rows, :])
    fg = seaborn.FacetGrid(data=df_sset, margin_titles=True, aspect=1,
            row='region', col='lam', hue='callset', sharey=True, hue_kws=dict(marker=__markers__))
    if vmc_pr is not None:
        lams = vmc_pr['lam'].cat.categories
        def helper(lamix):
            lam = lams[lamix]
            sel_rows = (vmc_pr['machine'] == 'Ada') \
                    & (vmc_pr['s2g'] == s2g) \
                    & (vmc_pr['sample'] == sample) \
                    & (vmc_pr['vartype'] == vartype) \
                    & (vmc_pr['lam'] == lam)
            def curveplotter(y='precision', linestyle='-', region='chr22'):
                df = vmc_pr.loc[sel_rows & (vmc_pr['region'] == region), :].copy()
                if region == 'chr22':
                    row = 2
                elif region == 'chr1_2':
                    row = 1
                fg.axes[row][lamix].plot(df['recall'], df[y],
                        color='black', linestyle=linestyle)
                return(None)
            curveplotter('precision', '-', region='chr22')
            curveplotter('precision_estim', ':', region='chr22')
            curveplotter('precision', '-', region='chr1_2')
            curveplotter('precision_estim', ':', region='chr1_2')
        [helper(ix) for ix in range(len(lams))]
    fg = fg.map(plt.plot, 'recall', 'precision')
    fg = fg.add_legend()
    return(fg)


def plotter2(df, hue='machine', sample='mix1'):
    '''
    Precision-recall plot; rows by s2g and columns by lambda
    '''
    seaborn.set()
    seaborn.set_context('notebook')
    sel_rows = (df['case_sample'] == sample)
    df_sset = df.loc[sel_rows, :]
    fg = seaborn.FacetGrid(data=df_sset, margin_titles=True, aspect=1,
            row='s2g', col='lam', hue=hue)
    fg = fg.map(plt.plot, 'recall', 'precision')
    #fg = fg.map(plt.plot, 'recall', 'precision_estim')
    fg = fg.add_legend()
    return(fg)


def plotter3(pr, vmc_pr=None, sample='mix1', region='autosomes', vartype='snp'):
    '''
    Precision-recall plot; rows by lambda and columns by s2g
    '''
    seaborn.set()
    seaborn.set_context('talk')
    sel_rows = (pr['sample'] == sample) & (pr['region'] == region) & (pr['vartype'] == vartype)
    df_sset = fix_names(pr.loc[sel_rows, :])
    fg = seaborn.FacetGrid(data=df_sset, margin_titles=True, aspect=1,
            row='lam', col='s2g', hue='callset', sharey=True,
            hue_kws=dict(marker=__markers__))
    fg = fg.map(plt.plot, 'recall', 'precision')
    fg = fg.add_legend()
    return(fg)


def plotter4(pr, vmc_pr=None, sample='mix1', lam=0.2, vartype='snp'):
    '''
    '''
    seaborn.set()
    seaborn.set_context('talk')
    # filter pr and vmc_pr
    pr = pr_astype(pr, False)
    # filter pr and vmc_pr
    pr = pr_astype(pr, False)
    vmc_pr = pr_astype(vmc_pr, True)
    sel_rows = (pr['sample'] == sample) & (pr['lam'] == lam) & (pr['vartype'] == vartype)
    df_sset = pr.loc[sel_rows, :]
    df_sset = pr_astype(df_sset, False)
    df_sset = fix_names(pr.loc[sel_rows, :])
    # all levels of regions in either pr or vmc_pr
    regions = vmc_pr['region'].cat.categories
    allregions = list(pr['region'].cat.categories)
    # create FacetGrid object without plotting
    fg = seaborn.FacetGrid(data=df_sset, margin_titles=True, aspect=1,
            row='region', col='s2g', hue='callset', sharey=True, hue_kws=dict(marker=__markers__))
    def helper(reg):
        '''
        Make all plots for a given region (a row of the plot matrix)
        '''
        # get index for reg
        regix = allregions.index(reg)
        sel_rows = (vmc_pr['machine'] == 'Ada') \
                & (vmc_pr['sample'] == sample) \
                & (vmc_pr['vartype'] == vartype) \
                & (vmc_pr['region'] == reg) \
                & (vmc_pr['lam'] == '0.2')
        def curveplotter(y='precision', linestyle='-', s2g=-3):
            # filter df
            df = vmc_pr.loc[sel_rows & (vmc_pr['s2g'] == s2g), :].copy()
            df = pr_astype(df, True)
            # determine column index based on s2g
            if s2g == -2:
                column = 2
            elif s2g == -3:
                column = 1
            elif s2g == -4:
                column = 0
            # plot on the axes object for row 'regix' and column 'column'
                column = 2
            elif s2g == -3:
                column = 1
            elif s2g == -4:
                column = 0
            # plot on the axes object for row 'regix' and column 'column'
            fg.axes[regix][column].plot(df['recall'], df[y],
                    color='black', linestyle=linestyle)
            return(None)
        curveplotter('precision', '-', s2g=-2)
        curveplotter('precision_estim', ':', s2g=-2)
        curveplotter('precision', '-', s2g=-3)
        curveplotter('precision_estim', ':', s2g=-3)
        curveplotter('precision', '-', s2g=-4)
        curveplotter('precision_estim', ':', s2g=-4)
        return(None)
    r = [helper(x) for x in regions]
    fg = fg.map(plt.plot, 'recall', 'precision')
    fg = fg.add_legend()
    return(fg)


def plotter5(pr, s2g=-3, region='autosomes', vartype='snp', onepanel=False):
    '''
    Precision-recall plot; hue by callset, rows by lambda and columns by s2g

    Parameters:
    pr: a precision recall data frame
    s2g: '-2', '-3' or '-4'
    region: chr22 or autosomes
    vartype: snp or indel
    onepanel: whether to draw only a single panel

    Returns:
    a FacetGrid plot object
    '''
    seaborn.set()
    size='notebook'
    row = 'lam'
    sel_rows = (pr['s2g'] == s2g) & (pr['region'] == region) & (pr['vartype']
            == vartype) & (pr['control_sample'] != 'no_ctr')
    if onepanel:
        row = None
        sel_rows = (pr['control_sample'] == 'mix3') & (pr['lam'] == 0.2) & sel_rows
        size='paper'
    seaborn.set_context(size)
    df_sset = pr.loc[sel_rows, :]
    fg = seaborn.FacetGrid(col='control_sample', aspect=1,
            row=row, hue='callset', data=df_sset,
            margin_titles=True, hue_kws=dict(marker=__markers__))
    fg = fg.map(plt.plot, 'recall', 'precision')
    fg = fg.add_legend()
    return(fg)


def plotter6(pr, region='autosomes', vartype='snp', explanvar='control_sample'):
    '''
    Precision-recall plot; hue by some explanatory variable, columns by callset

    Parameters:
    pr: a precision recall data frame
    region: chr22 or autosomes
    vartype: snp or indel
    explanvar: the explanatory variable whose effect we study; either lam, control_sample or s2g

    Returns:
    a FacetGrid plot object

    Details:
    Columns are wrapped.  strelka2Germline is excluded
    '''
    seaborn.set()
    seaborn.set_context('notebook')
    sel_rows =  (pr['region'] == region) & (pr['vartype'] == vartype) & \
            (pr['control_sample'] != 'no_ctr') & (pr['callset'] != 'strelka2Germline')
    if explanvar == 'control_sample':
        lam=0.2
        s2g=-3
        sel_rows = sel_rows & (pr['s2g'] == s2g) & (pr['lam'] == lam)
        marker = ['$1$', '$2$', '$3$']
    if explanvar == 'lam':
        control_sample='mix3'
        s2g=-3
        sel_rows = sel_rows & (pr['s2g'] == s2g) & (pr['control_sample'] == control_sample)
        # r for rapidly, s for slowly decaying exponential
        marker = ['$s$', '$r$', '$0$']
    if explanvar == 's2g':
        control_sample='mix3'
        lam=0.2
        sel_rows = sel_rows & (pr['lam'] == lam) & (pr['control_sample'] == control_sample)
        marker = ['$4$', '$3$', '$2$']
    df_sset = pr.copy().loc[sel_rows, :]
    df_sset['callset'] = df_sset['callset'].cat.remove_unused_categories()
    fg = seaborn.FacetGrid(col='callset', aspect=1,
            col_wrap=3, hue=explanvar, data=df_sset,
            hue_kws=dict(marker=marker))
    fg = fg.map(plt.plot, 'recall', 'precision', linestyle='')
    fg = fg.add_legend()
    return(fg)


def plotter7(df, otherdata=False):
    dt = df.copy()
    if not otherdata:
        dt = df.loc[df['callset'].isin(__callsets__), :].copy()
        dt['callset'] = dt['callset'].cat.set_categories(__callsets__)
    seaborn.set()
    seaborn.set_context('paper')
    g = seaborn.FacetGrid(hue='callset', data=dt, aspect=1,
            hue_kws=dict(marker=__markers__ + ['H', 'h', '+']))
    g = g.map(plt.plot, 'recall', 'precision', marker = 'o')
    g.add_legend()
    return(g)


def vmc_read_svmprob(vmcVCF):
    '''
    Read SVMPROB and other fields from a VCF produced by VMC.

    Unfortunately for some reason "bcftools normalize -m-" yields an error and
    truncates the VCF therefore multiallelic records cannot be split into
    multiple biallelic records.  The awkward workaround is to remove all but
    the first variant from the record with a sed command.

    Parameter:
    vmcVCF: path to the VCF

    Returns:
    a DataFrame whose columns correspond to VCF fields
    '''
    args0 = ['bcftools', 'query', '-f',
            '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVMPROB\n', vmcVCF]
    args1 = ['sed', 's/,\S\+//g']
    p0 = subprocess.Popen(args0, stdout=subprocess.PIPE)
    p1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stdin=p0.stdout)
    cnames = ['chrom', 'pos', 'ref', 'alt', 'svmprob']
    df = pd.read_csv(p1.stdout, sep='\t', header=None, names=cnames)
    return(df)


def nrecords_in_vcf(vcf, PASS=False):
    '''
    Count records in VCF

    Parameters:
    vcf: the path to the VCF

    Returns:
    the number of records
    '''
    if PASS:
        args0 = ['bcftools', 'view', '-f .,PASS', '--threads', __addthreads__, '-H', vcf]
    else:
        args0 = ['bcftools', 'view', '--threads', __addthreads__, '-H', vcf]
    args1 = ['wc', '-l']
    proc0 = subprocess.Popen(args0, shell=False, stdout=subprocess.PIPE)
    proc1 = subprocess.Popen(args1, shell=False, stdout=subprocess.PIPE, stdin=proc0.stdout)
    proc0.stdout.close()
    nrec = proc1.communicate()[0]
    nrec = int(nrec) # turn bytesliteral (e.g. b'5226\n') to integer
    return(nrec)


def vmc_precrecall(csetVCF, tsetVCF):
    '''
    Sort callset for svmprob and get precision-recall based on a truth set

    Parameters:
    csetVCF: path to the callset VCF
    tsetVCF: path to the truthset VCF

    Returns
    pr: a pandas DataFrame with precision and recall for sorted positions
    '''
    tempdir = tempfile.TemporaryDirectory()
    dname = tempdir.name
    falseposVCF = dname + os.path.sep + '0000.vcf.gz'
    falsenegVCF = dname + os.path.sep + '0001.vcf.gz'
    trueposVCF = dname + os.path.sep + '0002.vcf.gz'
    args0 = ['bcftools', 'isec', '-Oz', '-p', dname, csetVCF, tsetVCF]
    subprocess.run(args0)
    # read the positives into two DataFrames
    falsepos = vmc_read_svmprob(falseposVCF)
    truepos = vmc_read_svmprob(trueposVCF)
    falsepos['istrue'] = False
    truepos['istrue'] = True
    pr = pd.concat([falsepos, truepos])
    pr = pr.sort_values(by='svmprob', ascending=False)
    pr['TPsize'] = np.cumsum(pr['istrue'])
    pr['Psize'] = np.array(range(len(pr))) + 1
    pr['precision'] = pr['TPsize'] / pr['Psize'] 
    pr['recall'] = pr['TPsize'] / nrecords_in_vcf(tsetVCF)
    pr['precision_estim'] = np.cumsum(pr['svmprob']) / pr['Psize']
    return(pr)
