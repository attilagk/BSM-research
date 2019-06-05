import pandas as pd
import numpy as np
import tempfile
import shutil
import subprocess
import os
import re
import joint_gt_ceph as jgc
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
import functools
from scipy import stats

# number of additional threads to bcftools; greater than 0 slowed down
# ThinkStation so much that the Debian desktop became unresponsive
#__addthreads__ = str(os.cpu_count() - 3)
__addthreads__ = '0'

def make_ts_aaf(mix='mix1', vartype='snp', region='chr22', overwrite=True,
        tsdir='/home/attila/projects/bsm/results/2019-03-18-truth-sets/chr22/snp/truthset'):
    '''
    Create AAF-based truth set partitions given genotype-based partitions and
    count variants.

    Parameters:
    mix: mix1, mix2 or mix3
    vartype: snp or indel
    region: chr22 or autosomes
    overwrite: True or False depending on wheter to force overwriting existing files
    tsdir: the path to the main directory in the hierarchy

    Returns:
    a pandas data frame whose rows are truth subsets at various AAFs and whose
    "nvariants" column contains the number of variants
    '''
    def helper(aaf):
        genotypes = gt_of_aaf[mix][aaf]
        invcfs = [indir + os.sep + g + '.vcf.gz' for g in genotypes]
        outvcf = mixdir + os.sep + str(aaf) + '.vcf.gz'
        unsorted_outvcf = mixdir + os.sep + str(aaf) + '-unsorted.vcf.gz'
        args0 = ['bcftools', 'concat', '-a', '--threads', __addthreads__, '-o', unsorted_outvcf, '-Oz'] + invcfs
        args1 = ['bcftools', 'sort', '-o', outvcf, '-Oz', unsorted_outvcf]
        args2 = ['bcftools', 'index', '--threads', __addthreads__, '-t', outvcf]
        if overwrite or not os.path.isfile(outvcf):
            subprocess.run(args0)
            subprocess.run(args1)
            os.remove(unsorted_outvcf)
            subprocess.run(args2)
        # count records
        args3 = ['bcftools', 'view', '--threads', __addthreads__, '-H', outvcf]
        args4 = ['wc', '-l']
        proc1 = subprocess.Popen(args3, shell=False, stdout=subprocess.PIPE)
        proc2 = subprocess.Popen(args4, shell=False, stdout=subprocess.PIPE, stdin=proc1.stdout)
        proc1.stdout.close()
        nrec = proc2.communicate()[0]
        nrec = int(nrec) # turn bytesliteral (e.g. b'5226\n') to integer
        return(nrec)
    gt_of_aaf = jgc.get_gt_of_aaf()
    indir = tsdir + os.sep + 'genotypes'
    aafdir = tsdir + os.sep + 'aaf'
    if not os.path.isdir(aafdir):
        os.mkdir(aafdir)
    unfdir = aafdir + os.sep + 'unfiltered'
    if not os.path.isdir(unfdir):
        os.mkdir(unfdir)
    mixdir = unfdir + os.sep + mix
    if not os.path.isdir(mixdir):
        os.mkdir(mixdir)
    aafs = gt_of_aaf[mix].keys()
    nrec = [helper(f) for f in aafs]
    res = pd.DataFrame({'region': region, 'vartype': vartype, 'sample': mix, 'AAF': list(aafs),
        'count': nrec})
    res = res.astype({'region': 'category', 'vartype': 'category', 'sample': 'category'})
    return(res)

def make_ts_aaf_get_nvariants():
    '''
    Calls make_ts_aaf for all combination of regions, vartypes and samples

    Returns:
    a pandoc data frame with the number of variants at all those combinations and AAFs
    '''
    def tsaaf(mix, vartype, region):
        bdir = '/home/attila/projects/bsm/results/2019-03-18-truth-sets/'
        res = make_ts_aaf(mix=mix, vartype=vartype, region=region, tsdir=bdir + os.sep + region + os.sep + vartype + '/truthset')
        return(res)
    regions = ['chr22', 'autosomes']
    vartypes = ['snp', 'indel']
    samples = ['mix1', 'mix2', 'mix3']
    nvariants = pd.concat([tsaaf(mix=m, vartype=v, region=r) for r in regions
        for v in vartypes for m in samples])
    nvariants = nvariants.astype({'region': 'category', 'vartype': 'category', 'sample': 'category'})
    return(nvariants)

def evaluate_model(Y, xset, model, only_params):
    '''
    Evaluate a model at a set of x values (AAF values)

    Parameters:
    Y: total number of somatic variants
    xset: the set of x values (AAF values) at which somatic variants exist in the truth set
    model: L1, L2,...
    only_params: return only model parameters a and b

    Returns
    if only_params=False: the set of y values (number of modeled somatic variants)
    if only_params=True: the model parameters a (intercept) and b (slope)
    '''
    if model == 'L1':
        x_max = 50
    elif model == 'L2':
        x_max = 100
    elif model == 'L3':
        x_max = 50
    elif model == 'L4':
        x_max = 100
    xsubset = [z for z in xset if z <= x_max]
    n = len(xsubset)
    xsum = sum(xsubset)
    a = Y / (n - xsum / x_max)
    b = - a / x_max
    if model in ['L3', 'L4']:
        a = Y / n
        b = 0
    if only_params:
        return({'a': a, 'b': b})
    def L_function(x):
        y = a + b * x
        return(y)
    yset = [L_function(z) for z in xsubset] + [0] * (len(xset) - n)
    return(yset)


def evalmodel2df(nvariants, sample, vartype, region, model, Y, p_som2germ=None):
    '''
    Calls evaluate_model for a combination of Y (total n.o. variants), sample,
    vartype, region, model

    Parameters:
    nvariants: a data frame returned by make_ts_aaf_get_nvariants
    Y: total number of somatic variants
    sample: mix1, mix2 or mix3
    vartype: snp or indel
    region: chr22 or autosomes
    model: L1, L2,...
    p_som2germ: -log_10 of the ratio of the number of somatic to germline variants

    Returns:
    a pandas data frame with the modeled number of variants
    '''
    xset = nvariants.loc[(nvariants['region'] == region) &
            (nvariants['sample'] == sample) &
            (nvariants['vartype'] == vartype), 'AAF']
    y = evaluate_model(Y=Y, xset=xset, model=model, only_params=False)
    df = pd.DataFrame({'y': y, 'x': xset, 'model': model, 'Y': int(Y),
        'p_som2germ': p_som2germ, 'sample': sample, 'vartype': vartype, 'region': region})
    return(df)


def combine_regions_germ_vars(regions={'autosomes': 2929051733, 'chr22': 51304566},
        germ_vars = {'snp': 4e6, 'indel': 0.8e6}):
    '''
    Scales the number of germline variants according to the length of genomic
    regions

    Parameters:
    regions: a dictionary of lengths in bases for various regions
    germ_vars: the approximate expected number of germline snps and indels

    Returns:
    a pandoc data frame of approx number of germline variants at all combinations of region and variant type
    '''
    gk = germ_vars.keys()
    rk = regions.keys()
    res = {v: [int(germ_vars[v] * regions[r] / regions['autosomes']) for r in rk] for v in gk}
    res = pd.DataFrame(res, index=rk)
    return(res)


def evalmodel2df_all(nvariants, germ_vars=combine_regions_germ_vars()):
    '''
    Calls evalmodel2df for all combination of regions, vartypes samples,
    models and total n.o. variants

    Returns:
    a pandoc data frame with the number of modeled variants at all those combinations and AAFs
    '''
    regions = ['chr22', 'autosomes']
    vartypes = ['snp', 'indel']
    samples = ['mix1', 'mix2', 'mix3']
    models = ['L1', 'L2', 'L3', 'L4']
    p_som2germ = [2, 3, 4]
    l = [evalmodel2df(nvariants, sample=s, vartype=v, region=r, model=m,
        Y=germ_vars.at[r, v] * 10 ** (- s2g), p_som2germ=s2g)
            for s in samples for v in vartypes for r in regions for m in
            models for s2g in p_som2germ]
    return(pd.concat(l))

def get_taejeongs_aaf_sample(sample='S316', scale2pct=True):
    '''
    A single sample version of get_taejeongs_aaf; see details therein
    '''
    csv = '/big/data/bsm/Bae-2018-science/aan8690_TableS1/' + sample + '.csv'
    fr_cx = pd.read_csv(csv)['FR-CX']
    def helper(y):
        if isinstance(y, str) and re.search('somatic', y):
            def splitter(ix):
                return(np.int64(str(y).split(sep=':')[ix]))
            alt = splitter(-1)
            ref = splitter(-2)
            total = alt + ref
            return(alt / total)
        else:
            return(np.nan)
    aaf = [helper(y) for y in fr_cx]
    aaf = [y for y in aaf if not np.isnan(y)]
    if scale2pct:
        aaf = [100 * y for y in aaf]
    df = pd.DataFrame({'VAF': aaf, 'sample': sample})
    df = df.astype({'sample': 'category'})
    return(df)


def get_taejeongs_aaf(samples=['S316', 'S320'], scale2pct=True):
    '''
    Import Taejeong's VAF for all variants in a list of samples from his Science article

    See
    https://science.sciencemag.org/highwire/filestream/703017/field_highwire_adjunct_files/1/aan8690_TableS1.xlsx
    If multiple samples are given, a "pooled" sample is created from them.

    Parameters:
    samples: a list of samples (note that S275 fails to import)
    scale2pct: True if VAF should be in percent

    Returns:
    a pandas DataFrame with a column of VAF values and a sample column
    '''
    l = [get_taejeongs_aaf_sample(s, scale2pct=scale2pct) for s in samples]
    res = pd.concat(l)
    if len(samples) > 1:
        pooled = pd.DataFrame({'VAF': res['VAF'], 'sample': 'pooled'})
        res = pd.concat([res, pooled])
        res = res.astype({'sample': 'category'})
    #res = sum(l, [])
    return(res)


def lambda_hat(aaf):
    '''
    MLE of the rate parameter lambda from a sample of VAF

    Parameters:
    aaf: a list of VAF values

    Returns:
    MLE of lambda
    '''
    lhat = len(aaf) / sum(aaf)
    return(lhat)


def scaled_exponential1(lam, ntot):
    def expfun(aaf):
        aaf = np.array(aaf)
        y = lam * ntot * np.exp(- lam * aaf)
        return(y)
    return(expfun)


def scaled_exponential(lam, ntot, aafs):
    aafs = np.array(aafs) # ensure that aafs is a numpy array
    unscaled = lam * np.exp(- lam * aafs)
    scaled = ntot / sum(unscaled) * unscaled
    return(scaled)


def exp_model_df(nvariants, region, sample, vartype, log10s2g, lam):
    '''
    Create a histogram of AAF according to an exponential model

    Parameters:
    nvariants: the output of make_ts_aaf_get_nvariants, a pandas DataFrame
    region: 'autosomes' or 'chr22'
    sample: 'mix1', 'mix2', or 'mix3'
    vartype: 'snp' or 'indel'
    log10s2g: log base 10 of the odds of somatic : germline variants
    lam: rate lambda of the exponential

    Returns:
    a pandas DataFrame with AAF + count (i.e. histogram) and input parameters
    '''
    germ_vars = combine_regions_germ_vars()
    ntot = germ_vars[vartype][region] * 10 ** log10s2g
    ntot = np.int64(ntot)
    sel_rows = (nvariants['region'] == region) & (nvariants['sample'] ==
            sample) & (nvariants['vartype'] == vartype)
    aafs = nvariants.loc[sel_rows, 'AAF']
    count = scaled_exponential(lam, ntot, aafs)
    count = np.int64(count)
    d = {'AAF': aafs, 'count': count, 'region': region, 'sample': sample,
            'vartype': vartype, 'log10s2g': log10s2g, 'ntot': ntot, 'lambda': lam}
    df = pd.DataFrame(d)
    df = df.astype({'region': 'category', 'sample': 'category', 'vartype': 'category'})
    return(df)


def exp_model_df_concat(nvariants, log10s2gs=[-2, -3, -4], lambdas=[0.2, 0.04]):
    '''
    Concatenate into a DataFrame the iterator created by exp_model_iter

    Parameters:
    it: the iterator created by exp_model_iter
    log10s2gs: a list of multiple levels of log10s2g (see the log10s2g parameter of exp_model_df)
    lambdas: a list multiple levels of lambda (see the lam parameter of exp_model_df)

    Returns:
    a pandas DataFrame similar to the value of exp_model_df
    '''
    categcols = ['region', 'sample', 'vartype']
    l = [nvariants[c].cat.categories for c in categcols]
    l = l + [log10s2gs] + [lambdas]
    it = itertools.product(*l)
    df = pd.concat([exp_model_df(nvariants, *y) for y in it])
    df = df.astype({'lambda': 'category', 'log10s2g': 'category', 'region': 'category', 'sample': 'category', 'vartype': 'category'})
    return(df)


def exp_model_plot0(expm, log10s2g=-3, region='autosomes'):
    '''
    Plots a grid of histograms of AAF corresponding to a set of exponential models

    Parameters:
    expm: the value of exp_model_df_concat
    log10s2g: log base 10 of the odds of somatic : germline variants
    region: 'autosomes' or 'chr22'

    Returns:
    a seaborn FacetGrid object
    '''
    sel_rows = (expm['log10s2g'] == log10s2g) & (expm['region'] == region)
    df = expm.loc[sel_rows, :]
    g = sns.FacetGrid(df, row='sample', col='lambda', hue='vartype',
            sharey=False, aspect=2)
    g.map(plt.plot, 'AAF', 'count', marker='o', linestyle='dotted', markeredgecolor='white')
    g = g.add_legend()
    return(g)

def aaf_distplot(aafdf=get_taejeongs_aaf(), fit=stats.expon):
    '''
    Plot the distribution of VAF values 

    If the input aafdf contains multiple samples then a multi-plot grid is
    created.

    Parameters:
    aafdf: the pandas DataFrame output of get_taejeongs_aaf
    fit: None or a scipy.stats distribution

    Returns: a seaborn FacetGrid object
    '''
    sns.set()
    sns.set_context('talk')
    g =  sns.FacetGrid(aafdf, row='sample', aspect=2.5, height=4)
    g.map(sns.distplot, 'VAF', hist=True, rug=True, kde=False,
            fit=fit)
    if fit is None:
        pass
    else:
        pass
    return(g)


def aaf_distplot1(aafdf):
    samples = aafdf['sample'].cat.categories
    nsamples = len(samples)
    fig, ax = plt.subplots(nsamples, 1)
    def plot_sample(sample):
        ix = samples.searchsorted(sample)
        n, bins, patches = ax[ix].hist(aafdf[aafdf['sample'] == sample]['VAF'])
        res = pd.DataFrame({'count': list(n) + [list(n)[-1]], 'bins': bins,
            'sample': sample})
        ax[ix].set_title(sample)
        return(res)
    histo = [plot_sample(s) for s in samples]
    histo = pd.concat(histo)
    g = sns.FacetGrid(histo, row='sample')
    g.map(plt.step, 'bins', 'count', where='post')
    return(histo)



def downsample_aaf_vcf(ssize, invcfpath, outvcfpath, seed=19760415):
    '''
    Downsample one VCF

    Parameters:
    ssize: the number of calls in the sample (sample size)
    invcfpath: the path to the unfiltered input VCF
    outvcfpath: the path to the filtered output VCF
    seed: for sampling

    Returns:
    the regions (positions) of the sampled calls in a pandas DataFrame
    '''
    args0 = ['bcftools', 'view', '--threads', __addthreads__, '-H', invcfpath]
    args1 = ['cut', '-f1,2']
    proc0 = subprocess.Popen(args0, shell=False, stdout=subprocess.PIPE)
    proc1 = subprocess.Popen(args1, shell=False, stdout=subprocess.PIPE, stdin=proc0.stdout)
    regions = pd.read_csv(proc1.stdout, sep='\t', names=['CHROM', 'POS'], dtype={'CHROM': 'category', 'POS': np.int64})
    outs, errs = proc1.communicate()
    # downsampling regions
    sample = regions.sample(n=ssize, replace=False, axis=0, random_state=seed)
    sample = sample.sort_index()
    discarded_ix = set(regions.index) - set(sample.index)
    discarded = regions.loc[discarded_ix, :]
    discarded = discarded.sort_index()
    # output: regions file
    outdirpath = os.path.dirname(outvcfpath)
    outbname = os.path.basename(outvcfpath)
    discarded_vcfpath = outdirpath + os.sep + 'discarded-' + outbname
    if not os.path.exists(outdirpath):
        os.makedirs(outdirpath)
    if ssize == 0:
        shutil.copyfile(invcfpath, discarded_vcfpath)
        subprocess.run(['bcftools', 'index', '--tbi', discarded_vcfpath])
    else:
        regionspath = outvcfpath + '.regions'
        sample.to_csv(regionspath, sep='\t', header=False, index=False)
        # create reduced VCF
        args2 = ['bcftools', 'view', '--threads', __addthreads__, '-R', regionspath, '-Oz', invcfpath]
        args3 = ['bcftools', 'norm', '--rm-dup', 'both', '-Oz', '-o', outvcfpath]
        proc2 = subprocess.Popen(args2, shell=False, stdout=subprocess.PIPE)
        proc3 = subprocess.run(args3, shell=False, stdout=subprocess.PIPE, stdin=proc2.stdout)
        args4 = ['bcftools', 'index', '--tbi', outvcfpath]
        subprocess.run(args4)
        os.remove(regionspath)
        # create VCF from discarded records
        tempdir = tempfile.TemporaryDirectory()
        tempVCF = tempdir.name + os.path.sep + '0000.vcf.gz'
        args5 = ['bcftools', 'isec', '-C', '-Oz', '-p', tempdir.name, invcfpath, outvcfpath]
        subprocess.run(args5)
        shutil.move(tempVCF, discarded_vcfpath)
        shutil.move(tempVCF + '.tbi', discarded_vcfpath + '.tbi')
    return(sample)

def deduce_pathname(expm, topdir='/home/attila/projects/bsm/results/2019-03-18-truth-sets'):
    ix0 = expm.index[0]
    region = expm.at[ix0, 'region']
    vartype = expm.at[ix0, 'vartype']
    lam = expm.at[ix0, 'lambda']
    log10s2g = expm.at[ix0, 'log10s2g']
    sample = expm.at[ix0, 'sample']
    basedir = topdir + os.path.sep + region + os.path.sep + vartype + os.path.sep + 'truthset/aaf/'
    indir = basedir + 'unfiltered/' + sample + os.path.sep
    outdir = basedir + 'exp_model/lambda_' + str(lam) + os.path.sep + 'log10s2g_' + str(log10s2g) + os.path.sep + sample + os.path.sep
    d = {'indir': indir, 'outdir': outdir}
    return(d)


def downsample_all_aaf_vcfs(expm, topdir='/home/attila/projects/bsm/results/2019-03-18-truth-sets', seed=19760415):
    '''
    Downsample all VCFs referred to in expm

    Parameters:
    expm: a pandas DataFrame
    topdir: new directories and VCF-s will be created under it
    seed: for sampling

    Returns:
    the regions (positions) of the sampled calls in a dictionary of pandas
    DataFrames, where each key is an AAF
    '''
    d = deduce_pathname(expm=expm, topdir=topdir)
    indir = d['indir']
    outdir = d['outdir']
    outvcf = outdir + 'complete.vcf.gz'
    discarded_vcf = outdir + 'discarded-complete.vcf.gz'
    def helper(ix):
        '''
        Downsample a single VCF identified by ix of the pandas Data Frame expm
        '''
        aaf = expm.at[ix, 'AAF']
        ssize = expm.at[ix, 'count']
        invcfpath = indir + str(aaf) + '.vcf.gz'
        outvcfpath = outdir + str(aaf) + '.vcf.gz'
        discarded_vcfpath = outdir + 'discarded-' + str(aaf) + '.vcf.gz'
        df = downsample_aaf_vcf(ssize, invcfpath, outvcfpath, seed=seed)
        outpath = {'outvcfpath': outvcfpath, 'discarded_vcfpath': discarded_vcfpath}
        return(outpath)
    vcflist = [helper(ix) for ix in expm.index]
    outvcflist = [y['outvcfpath'] for y in vcflist if
            os.path.isfile(y['outvcfpath'])]
    discarded_vcflist = [y['discarded_vcfpath'] for y in vcflist]
    concat_vcfs(outvcf=outvcf, invcfs=outvcflist)
    concat_vcfs(outvcf=discarded_vcf, invcfs=discarded_vcflist)
    return(outvcf)


def concat_vcfs(outvcf, invcfs):
    '''
    Concatenate a list of VCFs into outvcf then sort and index that

    Parameters:
    outvcf: path to the output VCF
    invcfs: list of paths to the input VCFs

    Returns: None
    '''
    if os.path.exists(outvcf):
        os.remove(outvcf)
    # concatenate input VCFs and sort
    args0 = ['bcftools', 'concat', '-a', '--threads', __addthreads__] + invcfs
    args1 = ['bcftools', 'sort', '-Oz', '-o', outvcf]
    proc0 = subprocess.Popen(args0, shell=False, stdout=subprocess.PIPE)
    proc1 = subprocess.run(args1, shell=False, stdout=subprocess.PIPE,
            stdin=proc0.stdout)
    # make index
    args2 = ['bcftools', 'index', '-f', '--threads', __addthreads__, '-t', outvcf]
    subprocess.run(args2)
    return(None)


def split_up_expm(expm):
    '''
    Splits up expm according to the combination of the levels of categorical
    columns

    Parameters:
    expm: the output of make_ts_aaf_get_nvariants or exp_model_df

    Returns:
    the content of expm split into a list of DataFrames
    '''
    categcols = [c for c in expm.columns if expm[c].dtype.name == 'category']
    def helper(t):
        zip_l = list(zip(categcols, t))
        bool_l = [expm[y[0]] == y[1] for y in zip_l]
        sel_rows = functools.reduce(lambda a, b: a & b, bool_l)
        #sel_rows = bool_accumulate(bool_l)
        df = expm.loc[sel_rows, :]
        return(df)
    l = list(itertools.product(*[expm[c].cat.categories for c in categcols]))
    res = [helper(y) for y in l]
    return(res)


def downsample_absolutely_all_vcfs(expm, topdir='/home/attila/projects/bsm/results/2019-03-18-truth-sets', seed=19760415):
    '''
    Create truthset VCFs for all combinations of parameters respecting the directory tree created by make_ts_aaf

    Parameters:
    expm: a pandas DataFrame
    topdir: new directories and VCFs will be created under it
    seed: for sampling

    Returns:
    the pathnames of output VCFs
    '''
    l = split_up_expm(expm)
    outvcfs = [downsample_all_aaf_vcfs(df, topdir=topdir, seed=seed) for df in l]
    return(outvcfs)


def prec_recall_absolutely_all_vcfs(expm, topdir='/home/attila/projects/bsm/results/2019-03-18-truth-sets', seed=19760415):
    '''
    TODO
    '''
    l = split_up_expm(expm)
    l = l[0:2] # for testing
    outdirs = [deduce_pathname(expm=y, topdir=topdir) for y in l]
    return(outdirs)
    #prec-recall-vcf -t $truthsetvcf $outdir/*.vcf.gz > $prdir/prec-recall.csv
    truthsetvcf
    args = ['prec-recall-vcf', '-t', outvcf]
    subprocess.run(args)
    pass


def bool_accumulate(bool_l):
    '''
    Currently this function has a semantic bug.
    '''
    end = bool_l.pop()
    front = bool_l.copy()
    if len(front) == 1:
        val = front[0] & end
        return(val)
    else:
        return(bool_accumulate(front))
