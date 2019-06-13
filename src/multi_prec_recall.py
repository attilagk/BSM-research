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

__callsetmaindir__ = '/home/attila/projects/bsm/results/calls/benchmark-mix1-mix3/'
__truthsetmaindir__ = '/home/attila/projects/bsm/results/2019-03-18-truth-sets/'
__outmaindir__ = '/home/attila/projects/bsm/results/2019-05-02-make-truth-sets/'
__expmsubdir__ = 'truthset/aaf/exp_model/lambda_'
#__addthreads__ = str(os.cpu_count() - 3)
__addthreads__ = '7'
__allthreads__ = str(int(__addthreads__) + 1)


def getVCFpaths(callsetbn=None, region='chr22', vartype='snp', lam='0.04',
        log10s2g='-2', sample='mix1', callsetdir=None):
    '''
    Create pathname for various input, output, intermediate VCF files.

    Parameters:
    callsetbn: None, 'Tnseq.vcf.gz' or ['lofreqSomatic.vcf.gz', 'Tnseq.vcf.gz']
    region: chr22 or autosomes
    vartype: snp or indel
    lam: '0.04' or '0.2'
    log10s2g: '-2', '-3' or '-4'
    sample: 'mix1', 'mix2' or 'mix3'

    Returns:
    
    a dictionary of pathnames

    The keys of the dictionary are as follows:

    callset: the original callset without filtering

    prepared_callset_dir: directory of VCFs filtered for region and vartype
    
    prepared_callset: filtered for region and vartype
    
    reduced_truthset: the truthset according to the exp_model with parameters
    lam, log10s2g, sample
    
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
    subdir1 = region + os.path.sep + vartype + os.path.sep
    subdir2 = lam + '/log10s2g_' + log10s2g + os.path.sep + sample + os.path.sep
    if vartype == 'snp':
        alt_vartype = 'snvs'
    elif vartype == 'indel':
        alt_vartype = 'indels'
    # directories
    truthsetdir = __truthsetmaindir__ + subdir1 + __expmsubdir__ + subdir2
    if callsetdir is None:
        callsetdir = __callsetmaindir__ + alt_vartype + os.path.sep
    prepared_callset_dir = __outmaindir__ + subdir1
    reduced_callset_dir = __outmaindir__ + subdir1 + __expmsubdir__ + subdir2
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


def prepare4prec_recall(region='chr22', vartype='snp'):
    '''
    Run prepare4prec-recall shell script on initial callset VCFs for a given
    region and variant type

    Parameter(s):
    region: chr22 or autosomes
    vartype: snp or indel

    Returns: a list of the pathname of output VCFs
    '''
    if vartype == 'snp':
        callers = ['lofreqSomatic', 'somaticSniper', 'strelka2Germline2s', 'strelka2Somatic', 'Tnseq']
    elif vartype == 'indel':
        callers = ['strelka2Germline2s', 'strelka2Somatic', 'Tnseq']
    callsetbn = [c + '.vcf.gz' for c in callers]
    VCFpaths = getVCFpaths(callsetbn=callsetbn, region=region, vartype=vartype)
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
    elif region == 'chr22':
        r_opt = ['-r22']
    args1 = ['prepare4prec-recall', '-p', __allthreads__, '-r22', '-v',
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
        log10s2g='-2', sample='mix1', overwrite=False):
    '''
    Reduces (discards nonvariants from) the prepared callsets for a given region, vartype and exp_model

    Parameters:
    callsetbn: None or a list of callset VCF basenames like ['Tnseq.vcf.gz',...]
    region: chr22 or autosomes
    vartype: snp or indel
    lam: '0.04' or '0.2'
    log10s2g: '-2', '-3' or '-4'
    sample: 'mix1', 'mix2' or 'mix3'
    overwrite: whether to overwrite existing reduced callsets

    Returns: a list of pathnames of the reduced callsets
    '''
    VCFpaths = getVCFpaths(region=region, vartype=vartype)
    if callsetbn is None:
        callsetbn = glob.glob(VCFpaths['prepared_callset_dir'] + '*.vcf.gz')
        callsetbn = [os.path.basename(y) for y in callsetbn]
        return(callsetbn)
    VCFpaths = getVCFpaths(callsetbn=callsetbn, region=region,
            vartype=vartype, lam=lam, log10s2g=log10s2g, sample=sample)
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
    #args = ['prec-recall-vcf', '-p', '1', '-t', truthset] + callsets
    proc1 = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    prcsv = pd.read_csv(proc1.stdout)
    return(prcsv)


def reduce_precrecall(region='chr22', vartype='snp', lam='0.04',
        log10s2g='-2', sample='mix1'):
    VCFpaths = reduce_prepared_callsets(region=region, vartype=vartype, lam=lam,
            log10s2g=log10s2g, sample=sample)
    #return(VCFpaths)
    truthset = VCFpaths['reduced_truthset']
    callsets = VCFpaths['reduced_callset']
    pr = prec_recall_one_truthset(truthset=truthset, callsets=callsets)
    pr['region'] = region
    pr['vartype'] = vartype
    pr['lam'] = lam
    pr['log10s2g'] = log10s2g
    pr['sample'] = sample
    pr = pr_astype(pr)
    return(pr)


def prepare_reduce_precrecall(region='chr22', vartype='snp'):
    '''
    Prepare and reduce callset and calculate precision and recall for a given
    region and variant type
    '''
    val = prepare4prec_recall(region=region, vartype=vartype)
    def process1exp_model(lam, log10s2g, sample):
        pr = reduce_precrecall(region=region, vartype=vartype, lam=lam,
                log10s2g=log10s2g, sample=sample)
        return(pr)
    lams = ['0.04', '0.2']
    log10s2gs = ['-2', '-3', '-4']
    samples = ['mix1', 'mix2', 'mix3']
    l = [process1exp_model(lam=l, log10s2g=g, sample=s) for l in lams for g in
            log10s2gs for s in samples]
    pr = pd.concat(l)
    pr = pr_astype(pr)
    return(pr)


def vmc_prepare_reduce_precrecall(csetVCF, region='chr22', vartype='snp'):
    callsetbn = [os.path.basename(csetVCF)]
    VCFpaths = getVCFpaths(callsetbn=callsetbn, region=region, vartype=vartype)
    outdir = VCFpaths['prepared_callset_dir']
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outVCF = do_prepare4prec_recall(csetVCF=csetVCF, outdir=outdir,
            region=region, vartype=vartype, normalize=False, PASS=False, overwrite=True)
    # model specific operations
    lam = '0.04'
    log10s2g = '-2'
    sample = 'mix1'
    VCFpaths = reduce_prepared_callsets(callsetbn=callsetbn, region=region, vartype=vartype, lam=lam,
            log10s2g=log10s2g, sample=sample, overwrite=True)
    return(VCFpaths)


def run_all():
    '''
    Prepare and reduce callset and calculate precision and recall for all
    regions and variant types
    '''
    regions = ['chr22', 'autosomes']
    vartypes = ['snp', 'indel']
    l = [prepare_reduce_precrecall(region=r, vartype=v) for r in regions for v
            in vartypes]
    pr = pd.concat(l)
    pr = pr_astype(pr)
    return(pr)


def pr_astype(pr):
    '''
    Set data types for a precision recall data frame

    Parameters:
    pr: a precision recall data frame

    Returns: the data frame with the same data but corrected data types
    '''
    keys = ['callset', 'region', 'vartype', 'lam', 'log10s2g', 'sample']
    d = {k: 'category' for k in keys}
    pr = pr.astype(d)
    return(pr)


def read_pr_csv(csvpath):
    '''
    Read precision recall data from a CSV into a data frame and set data types
    '''
    pr = pd.read_csv(csvpath)
    pr = pr_astype(pr)
    return(pr)


def plotter1(df, sample='mix1', log10s2g=-2, vartype='snp'):
    '''
    Precision-recall plot; rows by log10s2g and columns by lambda
    '''
    seaborn.set()
    sel_rows = (df['sample'] == sample) & (df['log10s2g'] == log10s2g) & (df['vartype'] == vartype)
    df_sset = df.loc[sel_rows, :]
    fg = seaborn.FacetGrid(data=df_sset,
            row='region', col='lam', hue='callset', sharey=True)
    fg = fg.map(plt.plot, 'recall', 'precision', marker='o')
    fg = fg.add_legend()
    return(fg)


def vmc_read_svmprob(vmcVCF):
    args0 = ['bcftools', 'query', '-f',
            '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVMPROB\n', vmcVCF]
    args1 = ['sed', 's/,\S\+//g']
    p0 = subprocess.Popen(args0, stdout=subprocess.PIPE)
    p1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stdin=p0.stdout)
    cnames = ['chrom', 'pos', 'ref', 'alt', 'svmprob']
    df = pd.read_csv(p1.stdout, sep='\t', header=None, names=cnames)
    return(df)


def nrecords_in_vcf(vcf):
    args0 = ['bcftools', 'view', '--threads', __addthreads__, '-H', vcf]
    args1 = ['wc', '-l']
    proc0 = subprocess.Popen(args0, shell=False, stdout=subprocess.PIPE)
    proc1 = subprocess.Popen(args1, shell=False, stdout=subprocess.PIPE, stdin=proc0.stdout)
    proc0.stdout.close()
    nrec = proc1.communicate()[0]
    nrec = int(nrec) # turn bytesliteral (e.g. b'5226\n') to integer
    return(nrec)


def vmc_precrecall(csetVCF, tsetVCF):
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
    pos = pd.concat([falsepos, truepos])
    pos = pos.sort_values(by='svmprob', ascending=False)
    pos['TPsize'] = np.cumsum(pos['istrue'])
    pos['Psize'] = np.array(range(len(pos))) + 1
    pos['precision'] = pos['TPsize'] / pos['Psize'] 
    pos['recall'] = pos['TPsize'] / nrecords_in_vcf(tsetVCF)
    return(pos)
